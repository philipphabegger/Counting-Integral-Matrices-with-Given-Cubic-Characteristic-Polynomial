#!/usr/bin/env python
# coding: utf-8

# # Counting Integral Matrices with Given Cubic Characteristic Polynomial

# *P. Habegger*
# 
# ---

# ## Setup of Section 3

# In[2]:


R.<A1,A2,A3,B1,B2,B3,C1,C2,C3,L1,L2,L3,a,s,t,d,u,v,x,y,a1,a2,a3,r,l1,l2,l3> = PolynomialRing(QQ)


# In[3]:


Q = B1*B2*C1 + B3*C1^2 - B1^2*C2 - B1*C1*C3
R2 = (B1*B2*C2 - B2^2*C1 - B3 *C1 *C2 + B1 *C2* C3 )*L1 + (B1 *C2 - B2* C1 )*L2 - C1* L3 - B2^3* C1 + B1 *B2^2 *C2 - 2*B2* B3* C1* C2 + B1 *B3* C2^2 + B1* B2* C2 *C3 - B3 *C1 *C2* C3 + B1* C2 *C3^2
R3 = (-B2* B3* C1 + B1 *B3* C2 - B3 *C1 *C3 + B1 *C3^2 )*L1 + (-B3* C1 + B1* C3 )*L2 + B1* L3 - B2^2 *B3 *C1 + B1* B2* B3 *C2 - B3^2* C1* C2 - B2 *B3* C1 *C3 + 2*B1 *B3* C2 *C3 - B3 *C1* C3^2 + B1* C3^3


# In[4]:


pretty_print("The polynomial Q: ",Q)
pretty_print("The polynomial R_2: ",R2)
pretty_print("The polynomial R_3: ",R3)


# In[5]:


R2==(B1*B2*C2-B2^2*C1  - B3*C1*C2 + B1*C2*C3)*L1\
    +(B1*C2-B2*C1)*L2 \
    -C1*L3 \
    -B2^3*C1 + B1*B2^2*C2 - 2*B2*B3*C1*C2 + B1*B3*C2^2 + B1*B2*C2*C3 - B3*C1*C2*C3 + B1*C2*C3^2                        


# In[6]:


R3==(-B2*B3*C1 + B1* B3* C2 - B3* C1*C3 + B1* C3^2)*L1\
    +(-B3*C1 + B1*C3)*L2\
    +B1*L3\
    -B2^2*B3 *C1  + B1*B2*B3*C2 - B3^2* C1* C2 - B2* B3* C1*C3 + 2*B1*B3*C2*C3 - B3* C1 *C3^2 + B1*C3^3


# Make sure that we get the correct characteristic polynomial.

# In[7]:


matrix(3,3,[-L1-B2-C3,R2/Q,R3/Q,B1,B2,B3,C1,C2,C3]).charpoly()


# ## Section 3, leading up to Lemma 8

# We establish displayed equation (17).

# In[8]:


A = matrix(3,3,[A1,A2,A3,B1,B2,B3,C1,C2,C3])
A.characteristic_polynomial().coefficients()[0]
(R2-Q*A2).substitute({L1:A.characteristic_polynomial().coefficients()[2],\
                      L2:A.characteristic_polynomial().coefficients()[1],\
                      L3:A.characteristic_polynomial().coefficients()[0]}) == 0,\
(R3-Q*A3).substitute({L1:A.characteristic_polynomial().coefficients()[2],\
                      L2:A.characteristic_polynomial().coefficients()[1],\
                      L3:A.characteristic_polynomial().coefficients()[0]}) == 0


# We express the matrix $A$ as in Lemma 8.

# In[9]:


A = matrix(3,3,[a1,a2,a3,r*u,a+s*u,t*u+d*x,r*v,s*v+d*y,a+t*v])
pretty_print(A)


# ## Section 4, Lemma 11

# We extract the linear coefficients of the characteristic polynomial of $A$ as expressed above. Recall that $d=0$ in this case.

# In[10]:


A.charpoly().coefficients()[1].substitute({d:0,a1:-l1-2*a-s*u-t*v}) == \
(-l1-2*a-s*u-t*v)*(s*u+t*v+2*a)+a*s*u+v*a*t+a^2-(a2*u+a3*v)*r


# In[11]:


tmp1 = l2-A.charpoly().coefficients()[1].substitute({d:0,a1:-l1-2*a-s*u-t*v})


# In[12]:


tmp2 = (s*u + t*v)^2 + (l1 + 3*a)*(s*u + t*v) + 3*a^2 + 2*a*l1 + l2 + (a2*u + a3*v)*r


# In[13]:


tmp1 == tmp2


# ## Section 5, Lemma 13

# Now for the matrix $\widetilde A$

# In[14]:


Atilde = matrix(3,3,[a1,r*a2,r*a3,u,a+s*u,t*u+d*x,v,s*v+d*y,a+t*v])
pretty_print(Atilde)


# Next we insert $\widetilde A$ into $R_2$ to get $r_2$.

# In[29]:


r2 = R2.substitute({B1:Atilde[1,0],B2:Atilde[1,1],B3:Atilde[1,2],C1:Atilde[2,0],C2:Atilde[2,1],C3:Atilde[2,2]})

print("Degree with respect to s:",r2.degree(s))
print("Degree with repsect to t:",r2.degree(t))
print("Degree with repsect to r:",r2.degree(r))
print("S^2T term of r_2: ",r2.coefficient({s:2,t:1}))
print("S^2 term of r_2: ",r2.coefficient({s:2,t:0}).factor())
print("ST term of r_2: ",r2.coefficient({s:1,t:1}).factor())
print("S term of r_2: ",r2.coefficient({s:1,t:0}).factor())
print("T term of r_2: ",r2.coefficient({s:0,t:1}).factor())
print("d^3 of r_2: ",r2.coefficient({s:0,t:0,d:3}).factor())
print("d^2 of r_2: ",r2.coefficient({s:0,t:0,d:2}).factor())
print("d^1 of r_2: ",r2.coefficient({s:0,t:0,d:1}).factor())
print("d^0 of r_2: ",r2.coefficient({s:0,t:0,d:0}).factor())


# Now let's divide by the non-constant terms by $xv^2-yu^2 = 1$.

# In[30]:


print("S^2 term of r_2: ",r2.coefficient({s:2,t:0})/(x*v^2-y*u^2))
print("ST term of r_2: ",r2.coefficient({s:1,t:1})/(x*v^2-y*u^2))
print("S term of r_2: ",r2.coefficient({s:1,t:0})/(x*v^2-y*u^2))
print("T term of r_2: ",r2.coefficient({s:0,t:1})/(x*v^2-y*u^2))


# In[31]:


r3 = R3.substitute({B1:Atilde[1,0],B2:Atilde[1,1],B3:Atilde[1,2],C1:Atilde[2,0],C2:Atilde[2,1],C3:Atilde[2,2]})

print("Degree of r_3 with respect to s: ",r3.degree(s))
print("Degree of r_3 with repsect to t: ",r3.degree(t))
print("Degree of r_3 with repsect to r: ",r3.degree(r))
print("ST^2 term of r_3: ",r3.coefficient({s:1,t:2}))
print("T^2 term of r_3: ",r3.coefficient({s:0,t:2}).factor())
print("ST term of r_3: ",r3.coefficient({s:1,t:1}).factor())
print("S term of r_3: ",r3.coefficient({s:1,t:0}).factor())
print("T term of r_3: ",r3.coefficient({s:0,t:1}).factor())
print("d^3 term of r_3: ",r3.coefficient({s:0,t:0,d:3}).factor())
print("d^2 term of r_3: ",r3.coefficient({s:0,t:0,d:2}).factor())
print("d^1 term of r_3: ",r3.coefficient({s:0,t:0,d:1}).factor())
print("d^0 term of r_3: ",r3.coefficient({s:0,t:0,d:0}).factor())


# Now let's divide by the non-constant terms by $xv^2-yu^2 = 1$.

# In[28]:


print("T^2 term of r_3: ",r3.coefficient({s:0,t:2})/(x*v^2-y*u^2))
print("ST term of r_3: ",r3.coefficient({s:1,t:1})/(x*v^2-y*u^2))
print("S term of r_3: ",r3.coefficient({s:1,t:0})/(x*v^2-y*u^2))
print("T term of r_3: ",r3.coefficient({s:0,t:1})/(x*v^2-y*u^2))

