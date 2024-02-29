#ifndef GZFILE_HH_
#define GZFILE_HH_

#include <cassert>
#include <string>
#include <zlib.h>
#include <tuple>

namespace {

class GZFile {
public:
	explicit GZFile(const char *path, const char *mode = "r")
		: _file{ ::gzopen(path, mode) }
	{}

	explicit GZFile(int fd, const char *mode = "r")
		: _file{ ::gzdopen(fd, mode) }
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
	
	bool is_error() const
	{
		int err;
		return is_open() && ::gzerror(_file, &err) == nullptr;
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

	bool put(char ch)
	{
		assert(is_open());
		return ::gzputc(_file, ch) != -1;
	}

	bool putline(const std::string& line)
	{
		for (char ch : line)
			if (!put(ch))
				return false;

		return put('\n');
	}

private:
	GZFile(const GZFile&) = delete;
	GZFile& operator=(const GZFile&) = delete;

	gzFile _file;

}; // GZFile

} // (anonymous)

#endif // GZFILE_HH_
