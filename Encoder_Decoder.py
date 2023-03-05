'''
Nik Jensen
ECE 5660
--------------
Differential Endocder and Decoder for QPSK. This is used to correct 
phase ambiguity, where the carrier phase PLL can lock in place with 
incorrect phase offset. This could mean locking 90 degrees or more
off from the correct symbol. 
'''

def encode_qpsk_table(bits):
    print("Encoding given data")
    out_bits = 0b00
    #if first bit == 1 (i.e. 1011 or 1000 or 1100), flip all bits
    if (bits >= 8):
        print("x = ",bin(bits))
        bits = (~bits & 0xF) #added 0xF due to ~bits flips to 2's complement
        print("~x = ",bin(bits))
    #if 0000 or 0110, 00
    if bits == 0b0000 or bits == 0b0110:
        out_bits = 0b00
    #if 0001 or 0100, 01
    if bits == 0b0001 or bits == 0b0100:
        out_bits = 0b01
    #if 0010 or 0111, 10
    if bits == 0b0010 or bits == 0b0111:
        out_bits = 0b10
    #if 0011 or 0101, 11
    if bits == 0b0011 or bits == 0b0101:
        out_bits = 0b11
    return out_bits

def decode_qpsk(bits):
    # print("Decoding given data")
    out_bits = 0b00
    #if first bit == 1 (i.e. 1011 or 1000 or 1100), flip all bits
    if (bits >= 8):
        print("x = ",bin(bits))
        bits = (~bits & 0xF) #added 0xF due to ~bits flips to 2's complement
        print("~x = ",bin(bits))
    #if 0000 or 0101, 00
    if bits == 0b0000 or bits == 0b0101:
        out_bits = 0b00
    #if 0001 or 0111, 01
    if bits == 0b0001 or bits == 0b0111:
        out_bits = 0b01
    #if 0010 or 0100, 10
    if bits == 0b0010 or bits == 0b0100:
        out_bits = 0b10
    #if 0011 or 0110, 11
    if bits == 0b0011 or bits == 0b0110:
        out_bits = 0b11
    return out_bits

def encode_qpsk(b1,b2,d1,d2):
    # print("encode")
    strb = str(b1) + str(b2) + str(d1) + str(d2)
    x = int(strb,2)
    est = bin(encode_qpsk_table(x)) #encoder table
    # print("Encoder ouptut = ",est)
    if(len(est[2:]) == 1):
        est.format(est, '#2b')
        strest = str(d1) + str(d2) + '0' + str(est[2:])
    else:
        strest = str(d1) + str(d2) + str(est[2:])
    return strest #returns string value of result

if __name__ == "__main__":
    print("Starting encoder/decoder process...\n")
    binary_input = input("input qpsk input (b_2k, b_2k+1, d_2k-2, d_2k-1): ")
    b1 = binary_input[0] #b_2k
    b2 = binary_input[1] #b_2k+1
    d1 = binary_input[2] #d_2k-2
    d2 = binary_input[3] #d_2k-1
    print("\nbits out should match 0b" + str(b1) + str(b2))
    strest = encode_qpsk(b1,b2,d1,d2)
    print("strest = ",strest)
    est_ = int(strest,2) #convert string to number for decoder
    print("est_ = ", est_)
    dec = decode_qpsk(est_) #decoder
    print("results = ", bin(dec))

    