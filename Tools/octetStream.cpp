#include <string>
#include <sodium.h>

#include "Tools/octetStream.h"
#include "Networking/Player.h"

/**
 * @brief Construct a new octet Stream from a string
 * 
 * @param other other string
 */
octetStream::octetStream(const string& other)
{
    mxlen = other.size();
    len = mxlen;
    ptr = 0;
    data = new octet[mxlen];
    avx_memcpy(data, (const octet*)other.data(), len*sizeof(octet));
}

/**
 * @brief Initial allocation
 * 
 * @param maxlen 
 */
octetStream::octetStream(size_t maxlen)
{
    mxlen = maxlen;
    len = 0;
    ptr = 0;
    data = new octet[mxlen];
}

octetStream::octetStream(size_t len, const octet* source):
    octetStream(len)
{
    append(source, len);
}

// Free the memory
void octetStream::clear()
{
    if(data){
        delete[] data;
    }
    data = 0;
    len = mxlen = ptr = 0;
}

// Assign with other octetStream
void octetStream::assign(const octetStream &os)
{
    if(os.len>mxlen){
        if(data){
            delete[] data;
        }
        mxlen = os.mxlen;
        data = new octet[mxlen];
    }
    len = os.len;
    memcpy(data, os.data, len*sizeof(octet));
    ptr = os.ptr;
}

/**
 * @brief Get string
 * 
 * @return string 
 */
string octetStream::str() const
{
    return string((char*)get_data(), get_length());
}

void octetStream::hash(octetStream& output) const
{
  assert(output.mxlen >= crypto_generichash_blake2b_BYTES_MIN);
  crypto_generichash(output.data, crypto_generichash_BYTES_MIN, data, len, NULL, 0);
  output.len=crypto_generichash_BYTES_MIN;
}

octetStream octetStream::hash() const
{
  octetStream h(crypto_generichash_BYTES_MIN);
  hash(h);
  return h;
}

bool octetStream::equals(const octetStream& a) const
{
  if (len!=a.len) { return false; }
  return memcmp(data, a.data, len) == 0;
}

void octetStream::append_random(size_t num)
{
  randombytes_buf(append(num), num);
}

void octetStream::concat(const octetStream& os)
{
  memcpy(append(os.len), os.data, os.len*sizeof(octet));
}

void octetStream::store_bytes(octet* x, const size_t l)
{
  encode_length(append(4), l, 4);
  memcpy(append(l), x, l*sizeof(octet));
}

void octetStream::get_bytes(octet* ans, size_t& length)
{
  auto rec_length = get_int(4);
  if (rec_length != length)
    throw runtime_error("unexpected length");
  memcpy(ans, consume(length), length * sizeof(octet));
}

void octetStream::store(int l)
{
  encode_length(append(4), l, 4);
}


void octetStream::get(int& l)
{
  l = get_int(4);
}

void octetStreams::reset(size_t num)
{
    resize(num);
    for(auto& o: *this){
        o.reset_write_head();
    }
}

// Free the memory of each octetStream
void octetStreams::clear()
{
    for(auto &o: *this){
        o.clear();
    }
}
