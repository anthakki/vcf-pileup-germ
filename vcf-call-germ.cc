
#include "gtbeta.h"
#include "vcfreader.hh"
#include <cmath>
#include <cstdio>
#include <cstdlib>
#include <list>
#include <map>
#include <unordered_set>

namespace {

static __attribute__((unused))
string_view
basename(string_view path)
{
	const char *p = &*path.begin();
	const char *q = &*path.end();

	while (p != q && q[-1] == '/')
		--q;
	for (const char *z = p; z != q; ++z)
		if (*z == '/')
			p = &z[1];

	return string_view{ p, q };
}

static __attribute__((unused))
std::pair<string_view, string_view>
splitext(string_view path)
{
	const char *p = path.begin();

	for (const char *z = p; z != path.end(); ++z)
		if (*z == '.')
			p = z;

	return std::make_pair( string_view{ path.begin(), p }, string_view{ p, path.end() } );
}

static
void
compute_GL_22(double *gl, const long *ad, const char *mut_mask, std::size_t samples, double we, double me, size_t at)
{
	std::size_t combs = samples > 0;
	for (std::size_t i = 0; i < samples; ++i)
		combs *= 2*(2+1)/2;

	std::vector<double> jl( combs );
	gtbeta_22_jl( jl.data(), reinterpret_cast<const unsigned long *>(ad), samples, we, me, at );

	for (std::size_t i = 0; i < samples; ++i)
		if (mut_mask[i])
			gtbeta_zap( jl.data(), samples, 2*(2+1)/2, i, 0 );

	gtbeta_22_gl( gl, jl.data(), samples );
}

static
int
compute_GL_22_and_bias(double *gl, const long *ad, const long *ad2, const char *mut_mask, std::size_t samples, double we, double me, size_t at)
{
	std::size_t combs = samples > 0;
	for (std::size_t i = 0; i < samples; ++i)
		combs *= 2*(2+1)/2;

	std::vector<double> jl( combs );
	gtbeta_22_jl( jl.data(), reinterpret_cast<const unsigned long *>(ad), samples, we, me, at );

	for (std::size_t i = 0; i < samples; ++i)
		if (mut_mask[i])
			gtbeta_zap( jl.data(), samples, 2*(2+1)/2, i, 0 );

	std::vector<double> jl2( combs*combs );
	gtbeta_22_jl( jl2.data(), reinterpret_cast<const unsigned long *>(ad2), 2*samples, we, me, at );

	for (std::size_t i = 0; i < samples; ++i)
		if (mut_mask[i])
			gtbeta_zap( jl2.data(), samples, 2*(2+1)/2 * 2*(2+1)/2, i, 0 );

	double p1 = -HUGE_VAL;
	for (std::size_t c = 0; c < combs; ++c)
		if (jl[c] > p1)
			p1 = jl[c];

	double p2 = -HUGE_VAL;
	for (std::size_t c = 0; c < combs; ++c)
		if (jl2[c] > p2)
			p2 = jl2[c];

	bool strand_bias = p2 > p1;

	gtbeta_22_gl( gl, jl.data(), samples );

	return strand_bias;
}

static
const char *
compute_GT_22(const double *gl)
{
	return gtbeta_22_gt(gl);
}

static
bool
gt_is_variant(const char *gt)
{
	// 0, 0/0, 0/0/0, ... => false , 0/1, 1/0, 1/1, 1|1, 1/1/1, 2/0/3 => true
	for (const char *p = gt; *p != '\0'; ++p)
		if ('0' <= *p && *p <= '9')
			if (*p > '0')
				return true;
			
	return false;
}

} // (anonymous)

int
main(int argc, char **argv)
{
	bool header = false;
	double me = 0.10;
	double we = -1.;
	bool check_bias = false;
	bool known_somatic = false;
	std::map<std::string, bool> germ_samples;
	bool germ_default = false;
	std::size_t approx_taper = std::size_t(-1);

	array_view<const char *const> args{ &argv[1], &argv[argc] };

	for (; args.end() - args.begin() >= 1 && ( *args.begin() )[0] == '-'; args = { args.begin() + 1, args.end() })
		     if (string_view( *args.begin() ).str() == "-h")
			header = true;
		else if (string_view( *args.begin() ).str() == "-e")
		{
			if (parse_float( me, string_view(*( args.begin() + 1 )) ))
				args = { args.begin() + 1, args.end() };
			else
				goto usage;
		}
		else if (string_view( *args.begin() ).str() == "-ew")
		{
			if (parse_float( we, string_view(*( args.begin() + 1 )) ))
				args = { args.begin() + 1, args.end() };
			else
				goto usage;
		}
		else if (string_view( *args.begin() ).str() == "-b")
			check_bias = true;
		else if (string_view( *args.begin() ).str() == "-Z")
			known_somatic = true;
		else if (string_view( *args.begin() ).str() == "-s")
		{
			std::vector<string_view> samples;

			tokenize(samples, string_view(*( args.begin() + 1 )), ',');
			for (const auto& sample : samples)
				germ_samples.insert(germ_samples.end(), std::make_pair(sample.str(), true));
			germ_default = false;

			args = { args.begin() + 1, args.end() };
		}
		else if (string_view( *args.begin() ).str() == "-S")
		{
			std::vector<string_view> samples;

			tokenize(samples, string_view(*( args.begin() + 1 )), ',');
			for (const auto& sample : samples)
				germ_samples.insert(germ_samples.end(), std::make_pair(sample.str(), false));
			germ_default = true;

			args = { args.begin() + 1, args.end() };
		}
		else if (string_view( *args.begin() ).str() == "-at")
		{
			if (parse_unsigned( approx_taper, string_view(*( args.begin() + 1 )) ))
				args = { args.begin() + 1, args.end() };
			else
				goto usage;
		}
		else
			goto usage;

	if (!(we > 0.))
		we = 1. / me;

	if (args.end() - args.begin() != 1)
	{
usage:
		std::fprintf(stderr, "Usage: %s [-h] [-e em] [-ew ew] [-b] [-Z] [-s a,...] [-S a,...] [-at at] input.vcf" "\n", argv[0]);
		return EXIT_FAILURE;
	}

	const char *vcf_fn = *args.begin();

	GZFile vcf_fp{ vcf_fn, "r" };
	if (!vcf_fp.is_open())
	{
		std::fprintf(stderr, "%s: %s: %s" "\n", argv[0], vcf_fn, "failed to open for read");
		return EXIT_FAILURE;
	}

	VCFReader vcf_reader{ vcf_fp };
	std::vector<bool> germ_sample_mask;

	{
		if (header)
		{
			std::fwrite( vcf_reader.meta().begin(), 1, vcf_reader.meta().end() - vcf_reader.meta().begin(), stdout );
			std::printf("##" "FORMAT=<ID=GT,Number=1,Type=String,Description=\"Genotype\">" "\n");
			std::printf("##" "FORMAT=<ID=GL,Number=G,Type=Float,Description=\"Genotype likelihoods\">" "\n");
			std::printf("##" "INFO=<ID=GERMLINE,Number=0,Type=Flag,Description=\"Call for a germline variant\">" "\n");

			if (check_bias)
				std::printf("##" "INFO=<ID=STRANDBIAS,Number=0,Type=Flag,Description=\"Call for a strand biased variant\">" "\n");

			// TODO: print program version/date
		}

		germ_sample_mask.clear();
		for (string_view sample : vcf_reader.sample_fields())
		{
			auto it = germ_samples.find(sample.str());
			germ_sample_mask.push_back( it != germ_samples.end() ? it->second : germ_default );
		}

		std::printf("%.*s\n", int( vcf_reader.line().end() - vcf_reader.line().begin() ), vcf_reader.line().begin());
	}

	while (vcf_reader.read( vcf_fp ))
	{
		VCFInfoString info{ vcf_reader.info() };
		VCFFormatString format{ vcf_reader.format() };

		std::vector<std::array<long, 2> > ad;
		std::vector<std::array<long, 2> > ad_fb;
		{
			for (string_view sample_data : vcf_reader.sample_fields())
			{
				format.set_data(sample_data);

				auto parse_counts = [](long v[2], string_view field) {
					std::vector<string_view> vals;
					tokenize(vals, field, ',');

					if (!parse_unsigned(v[0], vals.front()))
						v[0] = 0;

					{
						v[1] = 0;
					for (string_view val : array_view<string_view>{ &*vals.begin() + 1, &*vals.end() })
					{
						long v1;

						if (!parse_unsigned(v1, val))
						{
							v[1] = 0;
							break;
						}
						else 
							v[1] += v1;
					}
					}
				};

				ad.push_back( std::array<long, 2>{ 0, 0 } );
				parse_counts( ad.back().data(), format.get("AD")  );

				if (check_bias)
				{
					long adf[] = { 0, 0 };
					parse_counts( adf, format.get("ADF") );

					ad_fb.push_back( std::array<long, 2>{ adf[0], adf[1] } );
					ad_fb.push_back( std::array<long, 2>{ ad.back()[0] - adf[0], ad.back()[1] - adf[1] } );
				}
			}
		}

		std::vector<char> mut_mask;
		{ std::size_t i = 0;
		for (string_view field : vcf_reader.sample_fields())
		{
			std::ignore = field;
			mut_mask.push_back( known_somatic && !germ_sample_mask[i++] );
		}
		}

		std::vector<std::array<double, 2*(2+1)/2> > gl{ ad.size() };
		bool is_bias = false;
		if (!check_bias)
			compute_GL_22( &gl[0][0], &ad[0][0], mut_mask.data(), ad.size(), we, me, approx_taper );
		else
			is_bias = compute_GL_22_and_bias( &gl[0][0], &ad[0][0], &ad_fb[0][0], mut_mask.data(), ad.size(), we, me, approx_taper);

		/* for (std::size_t i = 0; i < gl.size(); ++i)
			std::fprintf(stderr, "sample #%zu GT=%s GL=%g,%g,%g AD=%lu,%lu\n", i+1, compute_GT_22( &gl[i][0] ), gl[i][0], gl[i][1], gl[i][2], ad[i][0], ad[i][1]); */

		bool is_germ = false;
		std::list<std::string> string_cache;

		{ std::size_t i = 0;
		for (string_view field : vcf_reader.sample_fields())
		{
			format.set_data(field);

			const char *gt = compute_GT_22( &gl[i][0] );

			if (germ_sample_mask[i] && gt_is_variant(gt))
				is_germ = true;

			string_cache.emplace_back( gt );
			format.set( "GT", string_cache.back() );

			char GL_str[1024];
			std::snprintf(GL_str, countof(GL_str), "%g,%g,%g", gl[i][0], gl[i][1], gl[i][2]);
			string_cache.emplace_back( GL_str );
			format.set( "GL", string_cache.back() );

			string_cache.emplace_back( format.data_str() );
			vcf_reader.set_sample_field( i++, string_cache.back() );
		}
		}

		string_cache.emplace_back( format.head_str() );
		vcf_reader.set_format( string_cache.back() );

		{
			if (is_germ)
				info.set( "GERMLINE", "" );
			else
				info.unset( "GERMLINE" );

			string_cache.emplace_back( info.str() );
			vcf_reader.set_info( string_cache.back() );
		}

		if (check_bias)
		{
			if (is_bias)
				info.set( "STRANDBIAS", "" );
			else
				info.unset( "STRANDBIAS" );

			string_cache.emplace_back( info.str() );
			vcf_reader.set_info( string_cache.back() );
		}

		{
		const char *sep = "";
		for (string_view field : vcf_reader.fields())
		{
			std::printf("%s%.*s", sep, int( field.end() - field.begin() ), field.begin());
			sep = "\t";
		}
		std::printf("\n");
		}
	}

	if (!vcf_fp.is_eof())
	{
		std::fprintf(stderr, "%s: %s: %s" "\n", argv[0], vcf_fn, vcf_reader.error());
		return EXIT_FAILURE;
	}

	return EXIT_SUCCESS;
}
