#ifndef BYTE_SWAP_HPP_
#define BYTE_SWAP_HPP_

template <typename T>
T ByteSwap(T t) {
  
  const int number_of_bytes = sizeof(T);

  union {
    T t;
    unsigned char bytes[number_of_bytes];
  } dat1, dat2;

  dat1.t = t;

  for (int id_byte = 0; id_byte < number_of_bytes / 2; ++id_byte)
    dat2.bytes[id_byte] = dat1.bytes[number_of_bytes - 1 - id_byte];

  return dat2.t;

}

#endif // BYTE_SWAP_HPP_
