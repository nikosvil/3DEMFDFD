import math
import numpy as np

def error(EtotalXreal, EtotalXimag, EtotalYreal, EtotalYimag, ETXreal, ETXimag,  ETYreal, ETYimag):

 abs_error_X_real = abs(EtotalXreal-ETXreal)
 abs_error_X_imag = abs(EtotalXimag-ETXimag)
 abs_error_Y_real = abs(EtotalYreal-ETYreal)
 abs_error_Y_imag = abs(EtotalYimag-ETYimag)

 file1=open('AbsErrorX_Real.txt','w') 
 file1.write(str(abs_error_X_real))
 file2=open('AbsErrorX_Imag.txt','w') 
 file2.write(str(abs_error_X_imag))
 file3=open('AbsErrorY_Real.txt','w') 
 file3.write(str(abs_error_Y_real))
 file4=open('AbsErrorY_Imag.txt','w') 
 file4.write(str(abs_error_Y_imag))

 rel_error_X_real=100*(abs(ETXreal-EtotalXreal)/abs(ETXreal))
 rel_error_X_imag=100*(abs(ETXimag-EtotalXimag)/abs(ETXimag))
 rel_error_Y_real=100*(abs(ETYreal-EtotalYreal)/abs(ETYreal))
 rel_error_Y_imag=100*(abs(ETYimag-EtotalYimag)/abs(ETYimag))

 #print(rel_error_X_real)

 file5=open('RelErrorX_Real.txt','w') 
 file5.write(str(rel_error_X_real))
 file6=open('RelErrorX_Imag.txt','w') 
 file6.write(str(rel_error_X_imag))
 file7=open('RelErrorY_Real.txt','w') 
 file7.write(str(rel_error_Y_real))
 file8=open('RelErrorY_Imag.txt','w') 
 file8.write(str(rel_error_Y_imag))

 
 file1.close()
 file2.close()
 file3.close()
 file4.close()
 file5.close()
 file6.close()
 file7.close()
 file8.close()
 return  


