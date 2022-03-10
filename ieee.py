# would take too long to implement
# and I have more things to do

def b2frac(bits):
  e = -1
  frac = 0
  for bit in bits:
    if bit not in {'0', '1'}:
      raise TypeError

    frac += int(bit)*2**e
    e -= 1
    
  return frac

def d2b(dec):
  int(dec)
  bit = ''
  
  while dec > 0:
    bit += str(dec % 2)
    dec = dec // 2
    
  bit = bit[::-1]
  return bit


def ff2b(dec):
  frac = '0.' + str(dec)
  frac = float(frac)

  bit = ''
  e = -1
  val = 0
  for i in range(23):
    aux = 2**e
    if val+aux > frac:
      bit += '0'
    else:
      bit += '1'
      val += aux
    
    e -= 1
    
  return bit
  
  
# # float 0.f
d = 6875
f = float(f'0.{d}')

# bits = '10011001100110011001101'
bits = ff2b(d)
conv = b2frac(bits)

ea = conv - f
er = abs(ea)/f 

print(f)
print(bits)
print(conv)
print(f'erro absoluto: {ea}')
print(f'erro relativo: {er}')
  
# aux = ''
# aux += str(6%2)
# aux += '1'
# aux += '1'
# print(aux)
# aux = aux[::-1]
# print(aux)