#ifndef ARRAY_VIEW_HH_
#define ARRAY_VIEW_HH_

#include <algorithm>

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
		{ return std::size_t( end() - begin() ); }

		int compare(const array_view& other) const
		{
			return std::lexicographical_compare( other.begin(), other.end(), begin(), end() ) -
				std::lexicographical_compare( begin(), end(), other.begin(), other.end() );
		}

	private:
		Type *_first, *_last;

	}; // array_view

} // (anonymous)

#endif // ARRAY_VIEW_HH_
