
#define countof(x) \
	(sizeof((x)) / sizeof(*(x)))

#include "fishertest.h"
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

	std::string str() const
		;

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
			return false;
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
		;

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

#include <cstdio>
#include <cstdlib>
#include <list>
#include <map>

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
