#ifndef STRING_VIEW_HH_
#define STRING_VIEW_HH_

#include "array_view.hh"
#include <cassert>
#include <cstring>
#include <string>

namespace {

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
			return begin();
		}

		std::string str() const
		{ return std::string{ begin(), end() }; }

	}; // string_view

} // (anonymous)

#endif // STRING_VIEW_HH_
