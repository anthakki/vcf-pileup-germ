#ifndef ORDEREDDICT_HH_
#define ORDEREDDICT_HH_

#include <cstddef>
#include <functional>
#include <utility>
#include <vector>

namespace {

	template<class Key, class Value, class Hash = std::hash<Key>, class KeyEqual = std::equal_to<Key> >
	class OrderedDict {
	public:
		typedef std::pair<const Key, Value> value_type;

		OrderedDict()
			: _data{}
			, _index{ }
			, _hash{}
			, _keyequal{}
		{
			_index.resize( _MIN_SLOTS, std::size_t(-1) );
		}

		~OrderedDict() = default;

		value_type *begin()
		{ return &_data[0]; }
		const value_type *begin() const
		{ return &_data[0]; }

		value_type *end()
		{ return &_data[_data.size()]; }
		const value_type *end() const
		{ return &_data[_data.size()]; }

		std::size_t size() const
		{ return _data.size(); }

		bool empty() const
		{ return _data.empty(); }

		template<class OtherKey>
		value_type *find(const OtherKey& key)
		{
			std::size_t slot = _slot(key);
			std::size_t i;

			for (i = slot; i < _index.size(); ++i)
				if (_index[i] == std::size_t(-1))
					return end();
				else if (_keyequal(_data[_index[i]].first, key))
					return &_data[_index[i]];
			for (i = 0; i < slot; ++i)
				if (_index[i] == std::size_t(-1))
					return end();
				else if (_keyequal(_data[_index[i]].first, key))
					return &_data[_index[i]];

			return end();
		}

		template<class OtherKey>
		const value_type *find(const OtherKey& key) const
		{ return const_cast<OrderedDict *>(this)->find(key); }

		template<class OtherKey, class OtherValue>
		std::pair<value_type *, bool> emplace(OtherKey&& key, OtherValue&& value)
		{
			std::size_t slot = _slot(key);
			std::size_t i;

			// TODO: Robin Hood would be nice but that needs memory..

			for (i = slot; i < _index.size(); ++i)
				if (_index[i] == std::size_t(-1))
					goto insert;
				else if (_keyequal(_data[_index[i]].first, key))
					return std::make_pair( &_data[_index[i]], false );
			for (i = 0; i < slot; ++i)
				if (_index[i] == std::size_t(-1))
					goto insert;
				else if (_keyequal(_data[_index[i]].first, key))
					return std::make_pair( &_data[_index[i]], false );

insert:
			_index[i] = _data.size();
			_data.emplace_back(key, value);

			if (!_is_valid_size(_data.size()))
				_rehash();

			return std::make_pair( &_data[_index[i]], true );
		}

	private:
		OrderedDict(const OrderedDict&) = delete;
		OrderedDict& operator=(const OrderedDict&) = delete;

		template<class OtherKey>
		std::size_t _slot(const OtherKey& key) const
		{
			return _hash(key) % _index.size();
		}

		bool _is_valid_size(std::size_t size) const
		{ 
			// NOTE: load factor of ~33%, expected lookup in 3.25 slots
			return size < _index.size() - _index.size() / 4;
		}

		// initial number of slots
		static const std::size_t _MIN_SLOTS = 24;

		void _rehash()
		{
			// NOTE: grow to 200%
			std::size_t new_size = _index.size() + _index.size() / 4;

			std::fill(_index.begin(), _index.end(), std::size_t(-1));
			_index.resize( new_size, std::size_t(-1) );

			for (const value_type *it = begin(); it != end(); ++it)
			{
				std::size_t slot = _slot(it->first);
				std::size_t i;

				for (i = slot; i < _index.size(); ++i)
					if (_index[i] == std::size_t(-1))
						goto insert;
				for (i = 0; i < slot; ++i)
					if (_index[i] == std::size_t(-1))
						goto insert;

insert:
				_index[i] = it - begin();
			}
		}

		std::vector<value_type> _data;
		std::vector<std::size_t> _index;
		Hash _hash;
		KeyEqual _keyequal;

	}; // OrderedDict

} // (anonymous)

#endif // ORDEREDDICT_HH_
