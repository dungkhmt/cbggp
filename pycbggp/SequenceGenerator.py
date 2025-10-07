import random as rd 

class SequenceGenerator:
 def __init__(self):
  return 
 def GenerateNonDecreasingSumToM(self,n,M,up,low=1): 
  #generate a sequence of n positive integers such that the sum is equal to M, sorted in increasing order
  a = [0 for _ in range(n)] 
  i = 0
  while i < n:
   t = M//(n-i)
   #print('i = ',i,'M = ',M,'low = ',low)
   if t < low:
    return None
   t = min(t,up)    
   a[i] = rd.randint(low,t)
   
   M = M - a[i]
   low = a[i]   
   i = i + 1
   #print('a[',i-1,' = ',a[i-1],'low = ',low,'M = ',M)
   
  return a 
  
 def GenerateIncreasingSumToM(self,n,M,up,low=1): 
  #generate a sequence of n positive integers such that the sum is equal to M, sorted in increasing order
  a = [0 for _ in range(n)] 
  i = 0
  while n > 0:
   t = (M-n*(n-1)//2)//n
   #print('i = ',i,'M = ',M,'low = ',low)
   if t < low:
    return None
   t = min(t,up)    
   a[i] = rd.randint(low,t)
   
   M = M - a[i]
   low = a[i]+1   
   i = i + 1
   n = n - 1
   #print('a[',i-1,' = ',a[i-1],'low = ',low,'M = ',M)
   
  return a 
  
SG = SequenceGenerator()
a = SG.GenerateNonDecreasingSumToM(10,100,20,3)
b = SG.GenerateIncreasingSumToM(10,1000,200,2)
print(a,b)  