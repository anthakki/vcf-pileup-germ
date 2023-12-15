
#define countof(x) \
	(sizeof((x)) / sizeof(*(x)))

#include "gtbeta.h"
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

template<class String>
static
std::vector<string_view>
tokenize(std::vector<string_view>& fields, const String& line, char sep)
{
	fields.clear();
	fields.emplace_back( &*line.begin(), &*line.begin() );

	for (const char& ch : line)
		if (ch == sep)
			fields.emplace_back( &(&ch)[1], &(&ch)[1] );
		else
			fields.back() = { fields.back().begin(), &(&ch)[1] };

	return fields;
}

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

	array_view<const string_view> fields() const
	{ return array_view<const string_view>{ &*_fields.begin(), &*_fields.end() }; }

	array_view<const string_view> sample_fields() const
	{ 
		if (_fields.size() < countof(_VCF_FIELDS))
			return array_view<const string_view>{};
		else
			return array_view<const string_view>{ &_fields[countof(_VCF_FIELDS)], &_fields[_fields.size()] };
	}

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

class VCFInfoString {
public:
	explicit VCFInfoString(string_view data = ".")
		: _fields{}
		, _index{}
	{
		if (data.compare(string_view{"."}) == 0)
			return;

		std::vector<string_view> fields;
		tokenize(fields, data, ';');

		for (string_view field : fields)
		{
			_fields.emplace_back( field, string_view{ field.end(), field.end() } );

			for (const char& ch : field)
				if (ch == '=')
				{
					_fields.back() = { string_view{ _fields.back().first.begin(), &ch },
						string_view{ &(&ch)[1], _fields.back().second.end()  } };
				}
		}

		for (const auto& field : _fields)
			_index.insert( _index.end(), std::make_pair( field.first, _index.size() ) );
	}

	~VCFInfoString() = default;

	std::string str() const
	{
		if (_fields.size() < 1)
			return ".";

		std::string str;
		const char *sep = "";

		for (const auto& field : _fields)
		{
			if ( field.second.end() - field.second.begin() > 0 )
				str += sep + field.first.str() + '=' + field.second.str();
			else
				str += sep + field.first.str();
			sep = ";";
		}

		return str;
	}

	array_view<const std::pair<string_view, string_view> > fields() const
	{ return array_view<const std::pair<string_view, string_view> >{ &*_fields.begin(), &*_fields.end() }; }

	string_view get(string_view key) const
	{
		auto it = _index.find(key);
		return it != _index.end() ? _fields[it->second].second : string_view{"."};
	}

	void set(string_view key, string_view value)
	{
		auto it = _index.find(key);
		if (it != _index.end())
			_fields[it->second].second = value;
		else
		{
			it = _index.insert( it, std::make_pair( key, _index.size()) );
			_fields.emplace_back( key, value );
		}
	}

private:
	VCFInfoString(const VCFInfoString&) = delete;
	VCFInfoString& operator=(const VCFInfoString&) = delete;

	struct _Less {
		bool operator()(string_view lhs, string_view rhs) const
		{ return lhs.compare(rhs) < 0; }
	};

	std::vector<std::pair<string_view, string_view> > _fields;
	std::map<string_view, std::size_t, _Less> _index;

}; // VCFInfoString

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

static
void
compute_GL_22(double *gl, const long *ad, std::size_t samples, double we, double me)
{
	std::size_t combs = 1;

	for (std::size_t i = 0; i < samples; ++i)
		combs *= 2*(2+1)/2;

	std::vector<double> jl( combs );
	gtbeta_22_jl( jl.data(), reinterpret_cast<const unsigned long *>(ad), samples, we, me );
	gtbeta_22_gl( gl, jl.data(), samples );
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

#include <cstdio>
#include <cstdlib>
#include <list>
#include <map>
#include <unordered_set>

int
main(int argc, char **argv)
{
	bool header = false;
	double me = 0.10;
	double we = -1.;
	std::map<std::string, bool> germ_samples;
	bool germ_default = false;

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
		else
			goto usage;

	if (!(we > 0.))
		we = 1. / me;

	if (args.end() - args.begin() != 1)
	{
usage:
		std::fprintf(stderr, "Usage: %s [-h] [-e em] [-ew ew] [-s a,...] [-S a,...] input.vcf" "\n", argv[0]);
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
			}
		}

		std::vector<std::array<double, 2*(2+1)/2> > gl{ ad.size() };
		compute_GL_22( &gl[0][0], &ad[0][0], ad.size(), we, me );

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

		if (is_germ)
		{
			info.set( "GERMLINE", "" );

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
