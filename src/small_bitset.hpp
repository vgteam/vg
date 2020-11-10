#ifndef VG_SMALL_BITSET_INCLUDED
#define VG_SMALL_BITSET_INCLUDED

#include <cassert>

#include <sdsl/bits.hpp>

namespace vg {

/**
 * A small bitset. We expect that the universe size is usually at most 64.
 */
class SmallBitset {

    public:
        typedef std::uint64_t value_type;

        SmallBitset() : universe_size(0), data({ static_cast<value_type>(0) }) {}

        explicit SmallBitset(size_t n) : universe_size(n) {
            if (this->small()) {
                this->data.value = 0;
            } else {
                size_t sz = this->data_size();
                this->data.pointer = new value_type[sz]();
            }
        }

        ~SmallBitset() {
            this->clear();
        }

        SmallBitset(const SmallBitset& another) {
            this->copy(another);
        }

        SmallBitset(SmallBitset&& another) {
            this->move(another);
        }

        SmallBitset& operator=(const SmallBitset& another) {
            if (&another != this) {
                this->clear();
                this->copy(another);
            }
            return *this;
        }

        SmallBitset& operator=(SmallBitset&& another) {
            if (&another != this) {
                this->clear();
                this->move(another);
            }
            return *this;
        }

        size_t size() const { return this->universe_size; }
        bool small() const { return (this->size() <= VALUE_BITS); }
        size_t data_size() const { return (this->size() + VALUE_BITS - 1) / VALUE_BITS; }

        size_t count() const {
            if (this->small()) {
                return sdsl::bits::cnt(this->data.value);
            } else {
                size_t result = 0;
                size_t sz = this->data_size();
                for (size_t i = 0; i < sz; i++) {
                    result += sdsl::bits::cnt(this->data.pointer[i]);
                }
                return result;
            }
        }

        void insert(size_t i) {
            if (this->small()) {
                this->data.value |= static_cast<value_type>(1) << i;
            }
            else {
                this->data.pointer[i >> VALUE_SHIFT] |= static_cast<value_type>(1) << (i & VALUE_MASK);
            }
        }

        bool contains(size_t i) const {
            if (this->small()) {
                return (this->data.value & (static_cast<value_type>(1) << i));
            } else {
                return (this->data.pointer[i >> VALUE_SHIFT] & (static_cast<value_type>(1) << (i & VALUE_MASK)));
            }
        }

        void operator|=(const SmallBitset& another) {
            assert(this->size() == another.size());
            if (this->small()) {
                this->data.value |= another.data.value;
            } else {
                size_t sz = this->data_size();
                for (size_t i = 0; i < sz; i++) {
                    this->data.pointer[i] |= another.data.pointer[i];
                }
            }
        }

        void operator&=(const SmallBitset& another) {
            assert(this->size() == another.size());
            if (this->small()) {
                this->data.value &= another.data.value;
            } else {
                size_t sz = this->data_size();
                for (size_t i = 0; i < sz; i++) {
                    this->data.pointer[i] &= another.data.pointer[i];
                }
            }
        }

    private:
        size_t universe_size;
        union {
            value_type value;
            value_type* pointer;
        } data;

        void clear() {
            if (!this->small()) {
                delete[] this->data.pointer;
            }
            this->universe_size = 0;
            this->data.value = 0;
        }

        void copy(const SmallBitset& another) {
            this->universe_size = another.universe_size;
            if (this->small()) {
                this->data.value = another.data.value;
            } else {
                size_t sz = this->data_size();
                this->data.pointer = new value_type[sz];
                for (size_t i = 0; i < sz; i++) {
                    this->data.pointer[i] = another.data.pointer[i];
                }
            }
        }

        void move(SmallBitset& another) {
            this->universe_size = another.universe_size;
            if (this->small()) {
                this->data.value = another.data.value;
            } else {
                this->data.pointer = another.data.pointer;
                another.universe_size = 0;
                another.data.value = 0;
            }
        }

        constexpr static size_t VALUE_SHIFT = 6;
        constexpr static size_t VALUE_BITS = static_cast<size_t>(1) << VALUE_SHIFT;
        constexpr static size_t VALUE_MASK = VALUE_BITS - 1;
};

}

#endif // VG_SMALL_BITSET_INCLUDED
