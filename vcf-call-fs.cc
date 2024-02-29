
#include "fishertest.h"
#include "gzfile.hh"
#include "string_view.hh"
#include "vcfreader.hh"
#include <cmath>
#include <cstdio>
#include <cstdlib>
#include <list>
#include <map>

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
double
compute_FS(unsigned long a, unsigned long b, unsigned long c, unsigned long d)
{
	unsigned long t[] = { a, b, c, d };
	double fs = -10. * log10(fishertest(t, 0));

	return fs > 0. ? fs : 0.;  // NB. convert -0. to 0.
}

static
double
compute_SOR(double a, double b, double c, double d)
{
	// https://gatk.broadinstitute.org/hc/en-us/articles/360046786792-StrandOddsRatio

	double la = log(a + 1.), lb = log(b + 1.), lc = log(c + 1.), ld = log(d + 1.);

	double log_R = la-lb - ( lc-ld );
	double log_SR = log_R + log1p(exp( -2*log_R ));

	double log_RefR = -fabs( la-lb );
	double log_AltR = -fabs( lc-ld );

	return log_SR + ( log_RefR - log_AltR );
}

} // (anonymous)

int
main(int argc, char **argv)
{
	bool header = false;

	array_view<const char *const> args{ &argv[1], &argv[argc] };

	for (; args.end() - args.begin() >= 1 && ( *args.begin() )[0] == '-'; args = { args.begin() + 1, args.end() })
		     if (string_view( *args.begin() ).str() == "-h")
			header = true;
		else
			goto usage;

	if (args.end() - args.begin() != 1)
	{
usage:
		std::fprintf(stderr, "Usage: %s [-h] input.vcf" "\n", argv[0]);
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

	{
		if (header)
		{
			std::fwrite( vcf_reader.meta().begin(), 1, vcf_reader.meta().end() - vcf_reader.meta().begin(), stdout );
			std::printf("##" "INFO=<ID=FS,Number=1,Type=Float,Description=\"Fisher strand bias test\">" "\n");
			std::printf("##" "INFO=<ID=SOR,Number=1,Type=Float,Description=\"Symmetric odds ratio strand bias test\">" "\n");

			// TODO: print program version/date
		}

		std::printf("%.*s\n", int( vcf_reader.line().end() - vcf_reader.line().begin() ), vcf_reader.line().begin());
	}

	while (vcf_reader.read( vcf_fp ))
	{
		VCFInfoString info{ vcf_reader.info() };
		VCFFormatString format{ vcf_reader.format() };

		long tab[2][2] = { { 0 } };
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

				long ad[2], adf[2];
				parse_counts( ad,  format.get("AD")  );
				parse_counts( adf, format.get("ADF") );

				tab[0][0] += adf[0];         // ref+
				tab[0][1] += ad[0] - adf[0]; // ref-
				tab[1][0] += adf[1];         // alt+
				tab[1][1] += ad[1] - adf[1]; // alt-
			}
		}

		double FS  = compute_FS(  tab[0][0], tab[0][1], tab[1][0], tab[1][1] );
		double SOR = compute_SOR( tab[0][0], tab[0][1], tab[1][0], tab[1][1] );

		std::string FS_str = std::to_string(FS);
		info.set( "FS", FS_str );
		std::string SOR_str = std::to_string(SOR);
		info.set( "SOR", SOR_str );
		std::string info_str = info.str();
		vcf_reader.set_info( info_str );

		const char *sep = "";
		for (string_view field : vcf_reader.fields())
		{
			std::printf("%s%.*s", sep, int( field.end() - field.begin() ), field.begin());
			sep = "\t";
		}

		// std::printf(" -> %zu,%zu,%zu,%zu -> FS=%g SOR=%g", tab[0][0], tab[0][1], tab[1][0], tab[1][1], FS, SOR);
		std::printf("\n");
	}

	if (!vcf_fp.is_eof())
	{
		std::fprintf(stderr, "%s: %s: %s" "\n", argv[0], vcf_fn, vcf_reader.error());
		return EXIT_FAILURE;
	}

	return EXIT_SUCCESS;
}
