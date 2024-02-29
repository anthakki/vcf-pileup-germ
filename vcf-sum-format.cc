
#define countof(x) \
	(sizeof((x)) / sizeof(*(x)))

#include <cassert>
#include <cmath>
#include <string>
#include <zlib.h>
#include <tuple>

namespace {

class GZFile {
public:
	explicit GZFile(const char *path, const char *mode = "r")
		: _file{ ::gzopen(path, mode) }
	{}

	~GZFile()
	{
		if (is_open())
			close();
	}

	bool is_open() const
	{ return _file != nullptr; }

	void close()
	{
		assert(is_open());

		std::ignore = ::gzclose(_file);
		_file = nullptr;
	}

	bool is_eof() const
	{
		assert(is_open());
		return ::gzeof(_file) != 0;
	}

	int get()
	{
		assert(is_open());
		return gzgetc(_file);
	}

	bool getline(std::string& line)
	{
		line.clear();

		for (int ch; (ch = gzgetc(_file)) != -1;)
		{
			switch (ch)
			{
				case '\r':
					if ((ch = gzgetc(_file)) != -1)
						;
					else
						return is_eof();
					if (ch == '\n')
						;
					else
						return ::gzungetc(ch, _file) != -1;
					// FALLTHROUGH

				case '\n':
					return true;

				default:
					line.push_back(char(ch));
					break;
			}
		}

		return is_eof() && !line.empty();
	}

private:
	GZFile(const GZFile&) = delete;
	GZFile& operator=(const GZFile&) = delete;

	gzFile _file;

}; // GZFile

} // (anonymous)

#include <algorithm>
#include <cstring>
#include <string>

namespace {

template<class Type>
class array_view {
public:
	array_view()
		: array_view(nullptr, nullptr)
	{}

	array_view(const array_view&) = default;

	array_view(Type *first, Type *last)
		: _first{first}
		, _last{last}
	{}

	array_view& operator=(const array_view&) = default;

	~array_view() = default;

	Type *begin() const
	{ return _first; }

	Type *end() const
	{ return _last; }

	Type *data() const
	{ return begin(); }

	std::size_t size() const
	{ return end() - begin(); }

	int compare(const array_view& other) const
	{
		return std::lexicographical_compare( other.begin(), other.end(), begin(), end() ) -
			std::lexicographical_compare( begin(), end(), other.begin(), other.end() );
	}

private:
	Type *_first, *_last;

}; // array_view

class string_view : public array_view<const char> {
public:
	string_view() = default;

	string_view(const string_view&) = default;

	string_view(const char *first, const char *last)
		: array_view<const char>{ first, last }
	{}

	string_view(const std::string& string)
		: string_view{ string.data(), string.size() }
	{}

	string_view(const char *data)
		: string_view{ data, std::strlen(data) }
	{}

	string_view(const char *data, std::size_t size)
		: string_view{ &data[0], &data[size] }
	{}

	string_view& operator=(const string_view&) = default;

	~string_view() = default;

	const char *c_str() const
	{
		assert(*end() == '\0');
		return data();
	}

	std::string str() const
	{ return std::string{ begin(), end() }; }

}; // string_view

} // (anonymous)

#include <map>
#include <string>
#include <vector>

namespace {

template<class Unsigned, class String>
static
bool
parse_unsigned(Unsigned& value, const String& string)
{
	bool empty = true;

	value = 0;

	for (char ch : string)
	{
		if (!('0' <= ch && ch <= '9'))
			return false;

		Unsigned digit = ch - '0';
		if (!( value <= ( std::numeric_limits<Unsigned>::max() - digit ) / 10 ))
			return false;

		value = 10 * value + digit;
		empty = false;
	}

	return !empty;
}

template<class Signed, class String>
static
bool
parse_signed(Signed& value, const String& string)
{
	bool empty = true;

	value = 0;

	auto it = string.begin();
	for (; it != string.end(); ++it)
		if (*it == '-')
		{
			++it;
			break;
		}
		else
			return parse_unsigned(value, string);

	for (; it != string.end(); ++it)
	{
		if (!('0' <= *it && *it <= '9'))
			return false;

		Signed digit = *it - '0';
		if (!( value >= ( std::numeric_limits<Signed>::min() + digit ) / 10 ))
			return false;

		value = 10 * value - digit;
		empty = false;
	}

	return !empty;
}

template<class Float, class String>
static
bool
parse_float(Float& value, const String& string);

template<>
bool
parse_float<double, std::string>(double& value, const std::string& string)
{
	char *endp = nullptr;
	value = std::strtod( string.c_str(), &endp );
	return *endp == '\0';
}

template<>
bool
parse_float<double, string_view>(double& value, const string_view& string)
{ return parse_float(value, std::string{ string.begin(), string.end() }); }

template<class String, class Predicate>
static
std::vector<string_view>
tokenize(std::vector<string_view>& fields, const String& line, Predicate pred)
{
	fields.clear();
	fields.emplace_back( &*line.begin(), &*line.begin() );

	for (const char& ch : line)
		if (pred(ch))
			fields.emplace_back( &(&ch)[1], &(&ch)[1] );
		else
			fields.back() = { fields.back().begin(), &(&ch)[1] };

	return fields;
}

template<class String>
static
std::vector<string_view>
tokenize(std::vector<string_view>& fields, const String& line, char sep)
{ return tokenize(fields, line, [&](char ch) { return ch == sep; }); }

static const std::size_t _VCF_MIN_FIELDS = 8;
static const char *const _VCF_FIELDS[] = { "CHROM", "POS", "ID", "REF", "ALT", "QUAL", "FILTER", "INFO", "FORMAT" };

class VCFReader {
public:
	VCFReader()
		: _line{}
		, _fields{}
		, _meta{}
	{}

	explicit VCFReader(GZFile& file)
		: VCFReader{}
	{
		while (file.getline(_line))
			if (_line.size() > 0 && _line[0] == '#')
			{
				if (_line.size() > 1 && _line[1] == '#')
					_meta += _line + '\n';
				else
				{
					tokenize(_fields, string_view{ &*_line.begin() + 1, &*_line.end() }, '\t');

					if (!(_fields.size() >= _VCF_MIN_FIELDS &&
							std::equal( &_VCF_FIELDS[0], &_VCF_FIELDS[std::min( _fields.size(), countof(_VCF_FIELDS) )],
								_fields.begin(), [](string_view lhs, string_view rhs) { return lhs.compare(rhs) == 0; } )))
					{
						_throw_error("invalid VCF header");
						return;
					}

					break;
				}
			}
			else
			{
				_throw_error("missing VCF header");
				return;
			}
	}

	~VCFReader() = default;

	bool read(GZFile& file)
	{
		if (error())
			return false;

		while (file.getline(_line))
		{
			std::size_t width = _fields.size();

			tokenize(_fields, _line, '\t');
			if (_fields.size() < width)
				return _throw_error("missing VCF fields");
			if (_fields.size() > width)
				return _throw_error("excess VCF fields");

			return true;
		}

		if (!file.is_eof())
			return _throw_error("read error");

		return false;
	}

	const char *error() const
	{ return _fields.empty() ? _line.c_str() : nullptr; }

	string_view meta() const
	{ return _meta; }

	string_view line() const
	{ return _line; }

	string_view non_sample_line() const
	{
		auto fields = non_sample_fields();
		if (fields.size() > 0)
			return string_view{ fields.begin()->begin(), ( fields.end() - 1 )->end() };
		else
			return string_view{};
	}

	array_view<const string_view> fields() const
	{ return array_view<const string_view>{ &*_fields.begin(), &*_fields.end() }; }

	array_view<const string_view> sample_fields() const
	{ 
		if (_fields.size() < countof(_VCF_FIELDS))
			return array_view<const string_view>{ &_fields[_fields.size()], &_fields[_fields.size()] };
		else
			return array_view<const string_view>{ &_fields[countof(_VCF_FIELDS)], &_fields[_fields.size()] };
	}

	array_view<const string_view> non_sample_fields() const
	{ return array_view<const string_view>{ fields().begin(), sample_fields().begin() }; }

	void set_sample_field(std::size_t i, string_view field)
	{ _fields[ countof(_VCF_FIELDS) + i ] = field; }

	string_view chrom() const
	{ return _fields[0]; }

	string_view pos_str() const
	{ return _fields[1]; }

	long pos() const
	{
		unsigned long pos;

		if (parse_unsigned(pos, pos_str()))
			return long( pos - 1 );
		else
			return -1;
	}

	string_view id() const
	{ return _fields[2]; }

	string_view ref() const
	{ return _fields[3]; }

	string_view alt() const
	{ return _fields[4]; }

	string_view qual() const
	{ return _fields[5]; }

	string_view filter() const
	{ return _fields[6]; }

	string_view info() const
	{ return _fields[7]; }

	void set_info(string_view info)
	{ _fields[7] = info; }

	string_view format() const
	{
		if (_fields.size() < 9)
			return ".";

		return _fields[8];
	}

	void set_format(string_view format)
	{
		while (_fields.size() < 9)
			_fields.emplace_back(".");

		_fields[8] = format;
	}

protected:
	bool _throw_error(const char *error)
	{
		_line = error;
		_fields.clear();

		return false;
	}

private:
	VCFReader(const VCFReader&) = delete;
	VCFReader& operator=(const VCFReader&) = delete;

	std::string _line;
	std::vector<string_view> _fields;
	std::string _meta;

}; // VCFReader

class VCFFormatString {
public:
	explicit VCFFormatString(string_view head, string_view data = ".")
		: _fields{}
		, _index{}
	{
		if (head.compare(string_view{"."}) == 0)
			return;

		tokenize(_fields, head, ':');

		for (string_view field : _fields)
			_index.insert( _index.end(), std::make_pair( field, _index.size() ) );

		std::ignore = set_data(data);
	}

	~VCFFormatString() = default;

	std::string head_str()
	{
		std::vector<string_view> keys{ _fields.size(), "." };

		for (auto it = _index.begin(); it != _index.end(); ++it)
			keys[it->second] = it->first;

		if (_fields.size() < 1)
			return ".";

		std::string str;
		const char *sep = "";

		for (const auto& field : keys)
		{
			str += sep + field.str();
			sep = ":";
		}

		return str;
	}

	std::string data_str() const
	{
		if (_fields.size() < 1)
			return ".";

		std::string str;
		const char *sep = "";

		for (const auto& field : _fields)
		{
			str += sep + field.str();
			sep = ":";
		}

		return str;
	}

	bool set_data(string_view data)
	{
		if (data.compare(string_view{"."}) == 0)
		{
			std::fill( _fields.begin(), _fields.end(), "." );
			return true;
		}

		std::size_t width = _fields.size();

		tokenize(_fields, data, ':');
		if (_fields.size() < width)
		{
			while (_fields.size() < width)
				_fields.emplace_back(".");

			return false;
		}
		else if (_fields.size() > width)
			return false;

		return true;
	}

	string_view get(string_view key) const
	{
		auto it = _index.find(key);
		return it != _index.end() ? _fields[it->second] : string_view{"."};
	}

	void set(string_view key, string_view value)
	{
		auto it = _index.find(key);

		if (it != _index.end())
			_fields[it->second] = value;
		else
		{
			_fields.push_back(value);
			_index.insert( _index.end(), std::make_pair( key, _index.size() ) );
		}
	}

private:
	VCFFormatString(const VCFFormatString&) = delete;
	VCFFormatString& operator=(const VCFFormatString&) = delete;

	struct _Less {
		bool operator()(string_view lhs, string_view rhs) const
		{ return lhs.compare(rhs) < 0; }
	};

	std::vector<string_view> _fields;
	std::map<string_view, std::size_t, _Less> _index;

}; // VCFFormatString

static
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

static
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
