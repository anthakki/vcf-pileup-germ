#ifndef VCFREADER_HH_
#define VCFREADER_HH_

#define countof(x) \
	(sizeof((x)) / sizeof(*(x)))

#include "gzfile.hh"
#include "string_view.hh"
#include <map>
#include <string>
#include <vector>

namespace {

template<class Unsigned, class String>
static __attribute__((unused))
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
static __attribute__((unused))
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
static __attribute__((unused))
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
static __attribute__((unused))
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
static __attribute__((unused))
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

					return; // OK
				}
			}
			else
				break;

		_throw_error("missing VCF header");
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
			return string_view{ line().begin(), line().begin() };
	}

	string_view sample_line() const
	{
		auto fields = sample_fields();
		if (fields.size() > 0)
			return string_view{ fields.begin()->begin(), ( fields.end() - 1)->end() };
		else
			return string_view{ line().end(), line().end() };
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

	void unset(string_view key)
	{
		auto it = _index.find(key);
		if (it != _index.end())
		{
			_fields.erase( _fields.begin() + it->second );
			_index.erase(it);
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

} // (anonymous)

#endif // VCFREADER_HH_
