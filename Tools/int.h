#ifndef TOOLS_INT_H_
#define TOOLS_INT_H_

typedef unsigned char octet;
// Assumes word is a 64 bit value
#ifdef WIN32
  typedef unsigned __int64 word;
#else
  typedef unsigned long word;
#endif

/**
 * @brief positive modulo
 * 
 * @param i 
 * @param n 
 * @return int 
 */
inline int positive_modulo(int i, int n)
{
    return (i%n + n) % n;
}

inline void INT_TO_BYTES(octet *buff, int x)
{
    buff[0] = x&255;
    buff[1] = (x>>8)&255;
    buff[2] = (x>>16)&255;
    buff[3] = (x>>24)&255;
}

inline int BYTES_TO_INT(octet *buff)
{
    return buff[0] + 256*buff[1] + 65536*buff[2] + 16777216*buff[3];
}

#endif