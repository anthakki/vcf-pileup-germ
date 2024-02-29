#ifndef TEMPFILE_HH_
#define TEMPFILE_HH_

#include <string>
#include <unistd.h>

namespace {

	class TempFile {
	public:
		explicit TempFile(const char *path)
			: _path{path}
			, _fd{::mkstemp( &*_path.begin() )}
		{}

		~TempFile()
		{
			if (is_open())
				close();
		}

		bool is_open() const
		{ return _fd != -1; }

		bool close()
		{
			int result = ::close(_fd);
			_fd = -1;
			return result == 0;
		}

		int fd() const
		{ return _fd; }

		int release_fd()
		{
			int fd = _fd;
			_fd = -1;
			return fd;
		}

		const char *path() const
		{ return _path.c_str(); }

	private:
		TempFile(const TempFile&) = delete;
		TempFile& operator=(const TempFile&) = delete;

		std::string _path;
		int _fd;

	}; // TempFile

} // (anonymous)

#endif // TEMPFILE_HH_
