
#include "vcfreader.hh"

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

} // (anonymous)

#include <cstdio>
#include <cstdlib>
#include <list>

int
main(int argc, char **argv)
{
	bool header = false;
	std::vector<std::string> sample_names;

	array_view<const char *const> args{ &argv[1], &argv[argc] };

	for (; args.end() - args.begin() >= 1 && ( *args.begin() )[0] == '-'; args = { args.begin() + 1, args.end() })
		     if (string_view( *args.begin() ).compare(string_view("-h")) == 0)
			header = true;
		else if (string_view( *args.begin() ).compare(string_view("-S")) == 0)
		{
			std::vector<string_view> samples;

			tokenize(samples, string_view(*( args.begin() + 1 )), ',');
			for (const auto& sample : samples)
				sample_names.emplace_back(sample.str());

			args = { args.begin() + 1, args.end() };
		}
		else
			goto usage;

	if (args.end() - args.begin() < 1)
	{
usage:
		std::fprintf(stderr, "Usage: %s [-h] [-S a,...] input.vcf [...]" "\n", argv[0]);
		return EXIT_FAILURE;
	}

	std::list<GZFile> vcf_fps;
	std::list<VCFReader> vcf_readers;

	for (string_view vcf_fn : args) {
		vcf_fps.emplace_back( vcf_fn.c_str(), "r" );
		if (!vcf_fps.back().is_open())
		{
			std::fprintf(stderr, "%s: %s: %s" "\n", argv[0], vcf_fn.c_str(), "failed to open for read");
			return EXIT_FAILURE;
		}

		vcf_readers.emplace_back( vcf_fps.back() );
		if (vcf_readers.back().error())
		{
			std::fprintf(stderr, "%s: %s: %s" "\n", argv[0], vcf_fn.c_str(), vcf_readers.back().error());
			return EXIT_FAILURE;
		}
	}

	{
		for (const VCFReader& vcf_reader : vcf_readers)
			if (sample_names.empty())
				for (string_view sample : vcf_reader.sample_fields())
					sample_names.push_back(sample.str());
			else if (vcf_reader.sample_fields().size() != sample_names.size())
			{
				std::fprintf(stderr, "%s: %s: %s" "\n", argv[0], args.begin()[std::distance(vcf_readers.begin(), std::find_if(vcf_readers.begin(), vcf_readers.end(), [&](const VCFReader& lhs){ return &lhs == &vcf_reader; }))], "wrong number of samples");
				return EXIT_FAILURE;
			}

		if (header)
		{
			// TODO: print stuff from vcf_reader.meta()

			// TODO: print program version/date
		}

		{
			auto non_sample_line = string_view{ vcf_readers.front().line().begin(), vcf_readers.front().non_sample_line().end() };
			std::printf("%.*s", int(non_sample_line.size()), non_sample_line.data());

			for (const std::string& field : sample_names)
				std::printf("%c%.*s", '\t', int(field.size()), field.data());

			std::printf("\n");
		}
	}

	for (;;)
	{
		{
			auto vcf_fp_p = vcf_fps.begin();
			for (auto vcf_reader_p = vcf_readers.begin(); vcf_reader_p != vcf_readers.end(); ++vcf_reader_p, ++vcf_fp_p)
				if (!vcf_reader_p->read( *vcf_fp_p ))
				{
					if (vcf_fp_p == vcf_fps.begin())
					{
						for (++vcf_reader_p, ++vcf_fp_p; vcf_reader_p != vcf_readers.end(); ++vcf_reader_p, ++vcf_fp_p)
							if (vcf_reader_p->read( *vcf_fp_p ))
							{
								std::fprintf(stderr, "%s: %s: %s" "\n", argv[0], args.begin()[std::distance(vcf_readers.begin(), vcf_reader_p)], "excess data rows");
								return EXIT_FAILURE;
							}
					}
					else
							{
								std::fprintf(stderr, "%s: %s: %s" "\n", argv[0], args.begin()[std::distance(vcf_readers.begin(), vcf_reader_p)], "missing data rows");
								return EXIT_FAILURE;
							}

					goto done;
				}
		}

		for (const VCFReader& vcf_reader : vcf_readers)
			if (vcf_reader.non_sample_line().compare( vcf_readers.front().non_sample_line() ) != 0)
			{
				std::fprintf(stderr, "%s: %s: %s" "\n", argv[0], args.begin()[std::distance(vcf_readers.begin(), std::find_if(vcf_readers.begin(), vcf_readers.end(), [&](const VCFReader& lhs){ return &lhs == &vcf_reader; }))], "field data mismatch");
				return EXIT_FAILURE;
			}

		std::vector<std::string> new_sample_fields;

		for (const string_view& field : vcf_readers.front().sample_fields())
		{
			struct Data {
				enum { T_STRING, T_INT, T_DOUBLE } type;
				union {
					string_view v_string;
					double v_double;
				};
				string_view suffix;
			};

			std::vector<string_view> lhs_tokens;
			std::vector<Data> values;

			tokenize(lhs_tokens, field, [](char ch) { return ch == ':' || ch == ','; });
			for (std::size_t i = 0; i < lhs_tokens.size(); ++i)
			{
				string_view sfx{ lhs_tokens[i].end(), i < lhs_tokens.size() - 1 ? lhs_tokens[i+1].begin() : field.end() };
				long long v_int;
				double v_double;

				if (parse_signed(v_int, lhs_tokens[i]))
					values.push_back( Data{ .type = Data::T_INT, .v_double = double(v_int), .suffix = sfx });
				else if (parse_float(v_double, lhs_tokens[i]))
					values.push_back( Data{ .type = Data::T_DOUBLE, .v_double = v_double, .suffix = sfx });
				else
					values.push_back( Data{ .type = Data::T_STRING, .v_string = lhs_tokens[i], .suffix = sfx });
			}

			std::vector<string_view> rhs_tokens;

			for (const VCFReader& vcf_reader : vcf_readers)
				if (&vcf_reader != &vcf_readers.front())
				{
					if (tokenize(rhs_tokens, vcf_reader.sample_fields().begin()[ &field - vcf_readers.front().sample_fields().begin() ], [](char ch) { return ch == ':' || ch == ','; }).size() == lhs_tokens.size())
					{
						for (std::size_t i = 0; i < lhs_tokens.size(); ++i)
						{
							if (values[i].type == Data::T_STRING)
							{
								if (!(rhs_tokens[i].compare(lhs_tokens[i]) == 0))
								{
invalid_data:
									std::fprintf(stderr, "%s: %s: %s" "\n", argv[0], args.begin()[std::distance(vcf_readers.begin(), std::find_if(vcf_readers.begin(), vcf_readers.end(), [&](const VCFReader& lhs){ return &lhs == &vcf_reader; }))], "format data mismatch");
									return EXIT_FAILURE;
								}
							}
							else if (values[i].type == Data::T_INT)
							{
								long long v_int;
								double v_double;

								if (parse_signed(v_int, rhs_tokens[i]))
									values[i].v_double += v_int;
								else if (parse_float(v_double, rhs_tokens[i]))
									values[i] = Data{ .type = Data::T_DOUBLE, .v_double = values[i].v_double + v_double };
								else
									goto invalid_data;
							}
							else
							{
								double v_double;

								if (parse_float(v_double, rhs_tokens[i]))
									values[i].v_double += v_double;
								else
									goto invalid_data;
							}
						}
					}
					else
					{
						std::fprintf(stderr, "%s: %s: %s" "\n", argv[0], args.begin()[std::distance(vcf_readers.begin(), std::find_if(vcf_readers.begin(), vcf_readers.end(), [&](const VCFReader& lhs){ return &lhs == &vcf_reader; }))], "format structure mismatch");
						return EXIT_FAILURE;
					}
				}

			new_sample_fields.emplace_back();

			for (const Data& value : values)
			{
				if (value.type == Data::T_STRING)
					new_sample_fields.back().append( value.v_string.data(), value.v_string.size() );
				else if (value.type == Data::T_INT)
				{
					char buffer[1024];
					std::snprintf(buffer, sizeof(buffer), "%.0f", value.v_double);
					new_sample_fields.back().append( buffer );
				}
				else if (value.type == Data::T_DOUBLE)
				{
					char buffer[1024];
					std::snprintf(buffer, sizeof(buffer), "%g", value.v_double);
					new_sample_fields.back().append( buffer );
				}

				new_sample_fields.back().append( value.suffix.data(), value.suffix.size() );
			}
		}

		{
			std::printf("%.*s", int(vcf_readers.front().non_sample_line().size()), vcf_readers.front().non_sample_line().data());

			for (const std::string& field : new_sample_fields)
				std::printf("%c%.*s", '\t', int(field.size()), field.data());

			std::printf("\n");
		}
	}
done:
	;

	for (const GZFile& vcf_fp : vcf_fps)
		if (!vcf_fp.is_eof())
		{
			auto vcf_reader_p = vcf_readers.begin();
			std::advance(vcf_reader_p, std::distance(vcf_fps.begin(), std::find_if(vcf_fps.begin(), vcf_fps.end(), [&](const GZFile& lhs){ return &lhs == &vcf_fp; })));
			std::fprintf(stderr, "%s: %s: %s" "\n", argv[0], args.begin()[std::distance(vcf_fps.begin(), std::find_if(vcf_fps.begin(), vcf_fps.end(), [&](const GZFile& lhs){ return &lhs == &vcf_fp; }))], vcf_reader_p->error());
			return EXIT_FAILURE;
		}

	return EXIT_SUCCESS;
}
