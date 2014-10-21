

class HistogramSequence :
   def __init__(self,*arg,**karg) :
      """Creates an Histogram Sequence 
         HS=HistogramSequence() ; empty object
         HS=HistogramSequence(PickleFileName) ; load from a pickle file name
         HS=HistogramSequence(n_seq,minval,maxval,step) ; histogram for given set of minval maxval and so on
         HS=HistogramSequence(n_seq,base_array) ; histogram for given base (step will be computed from base_array)
         
         V0.9.2 - 2012 Aug 27 - 2012 Oct 05 - R.Gusteri, M.Maris
      """
      import numpy as np
      self._names_vectors=['H','quality','rows','rows_delta','stable_fraction','Ntot','Ninf','Nsup','time']
      self.kind='differential:number'
      self._empty_Iam()
      self._empty_Num()
      try :
         self.n_over=int(karg['n_over'])
      except :
         self.n_over=1
      try :
         self.rows=karg['rows']
      except :
         self.rows=None
      try :
         self.rows_delta=karg['rows_delta']
      except :
         self.rows_delta=None
      try :
         self.stable_fraction=karg['stable_fraction']
      except :
         self.stable_fraction=None
      try :
         self.heights=karg['heights']
      except :
         self.heights=None
      try :
         self.celltype=karg['celltype']
      except :
         self.celltype='uint'
      try :
         self.name=karg['name']
      except :
         self.name=None
      try :
         self.unit=karg['unit']
      except :
         self.unit=None
      try :
         self.row_name=karg['row_name']
      except :
         self.row_name=None
      try :
         self.row_unit=karg['row_unit']
      except :
         self.row_unit=None
      try :
         self.col_name=karg['col_name']
      except :
         self.col_name=None
      try :
         self.col_unit=karg['col_unit']
      except :
         self.col_unit=None
      if len(arg) == 1 :
         if type(arg[0]) == type(''):
            self.load(arg[0])
         return
      elif len(arg) == 2 :
         self.n_seq=arg[0]
         self.base=np.array(arg[1])
         self.n_steps=len(self.base)
         idx = np.where(abs(self.base) < np.inf)[0]
         self.minval=self.base[idx].min()
         self.maxval=self.base[idx].max()
         self.step=(self.maxval-self.minval)/float(len(idx))
      elif len(arg) == 4 :
         self.n_seq=arg[0]
         self.minval=arg[1]
         self.maxval=arg[2]
         self.step=arg[3]
         self.base = np.arange(self.minval-self.n_over*self.step-self.step,self.maxval+self.n_over*self.step+self.step+self.step,self.step)
         self.n_steps=len(self.base)
         self.base[0] = -np.inf
         self.base[-1] = +np.inf
         self.minval = self.minval - self.step*self.n_over
         self.maxval = self.maxval + self.step*self.n_over
      else :
         return 
      self.H=np.zeros([self.n_seq,self.n_steps],dtype=self.celltype)
      self.quality=np.zeros(self.n_seq,dtype='uint8')
   def _empty_Iam(self) :
      self.Iam={}
      self.Iam['differential']=True
      self.Iam['frequency']=False
   def _empty_Num(self) :
      self.Ninf=None
      self.Nsup=None
      self.Ntot=None
   def _recompute_base(self) :
      "used when base have to be recomputed, as an example after a col_slice"
      import numpy as np
      self.n_steps=len(self.base)
      idx = np.where(abs(self.base) != np.inf)[0]
      if len(idx) > 0 :
         self.minval = self.base[idx].min() - self.step*self.n_over
         self.maxval = self.base[idx].max() + self.step*self.n_over
   def _keys_vectors(self) :
      import numpy as np
      return np.array(self._names_vectors)
   def _make_copy_of_scalars(self) :
      import copy
      new = HistogramSequence()
      for k in self.__dict__.keys() :
         if (self._keys_vectors() == k).sum() == 0 :
            new.__dict__[k]=copy.deepcopy(self.__dict__[k])
      return new
   def _make_copy_of_all_but_H(self) :
      import copy
      new = HistogramSequence()
      for k in self.__dict__.keys() :
         if k != 'H' :
            new.__dict__[k]=copy.deepcopy(self.__dict__[k])
      return new
   def __len__(self) :
      return self.n_seq
   def version(self) :
      return ['1.0 - 28 aug 2012','1.1 - 10 set 2012 - ']
   def save(self,pickle_file) :
      import pickle
      if type(pickle_file)==type('') :
         try :
            pickle.dump(self.__dict__,open(pickle_file,'w'))
         except :
            return False
      else :
         try :
            pickle.dump(self.__dict__,pickle_file)
         except :
            return False
      return True
   def load(self,pickle_file) :
      import pickle
      if type(pickle_file)==type('') :
         try :
            self.__dict__=pickle.load(open(pickle_file,'r'))
         except :
            return False
      else :
         try :
            self.__dict__=pickle.load(pickle_file)
         except :
            return False
      if not self.__dict__.has_key('Iam') : self._empty_Iam()
      for k in ['Ntot','Ninf','Nsup'] : 
         if not self.__dict__.has_key(k) : 
            self._empty_Num()
      return True
   def copy(self) :
      import copy
      return copy.deepcopy(self)
   def valid_iseq(self,iseq):
      if iseq<-(self.n_seq-1) : 
         return False
      if iseq>(self.n_seq-1) :
         return False
      return True
   def fill(self,iseq,x,row_value=None,row_delta=None,stable_fraction=None,quality=None) :
      import numpy as np
      if self.valid_iseq(iseq)==False :
         raise Exception('iseq %d exceeds maximum lenght %d'%(iseq,len(self)))
      interval = np.array(((x-self.minval)/self.step).round(),dtype='int')+1
      idx=np.where(interval<=0)[0]
      if len(idx)>0 : interval[idx]=0
      idx=np.where(interval>=len(self.base)-1)[0]
      if len(idx)>0 : interval[idx]=len(self.base)-1
      for i in interval : self.H[iseq][i]+=1
      if self.rows != None and row_value != None :
         self.rows[iseq]=row_value
      if self.rows_delta != None and row_delta != None :
         self.rows_delta[iseq]=row_delta
      if self.stable_fraction != None and stable_fraction != None :
         self.stable_fraction[iseq]=stable_fraction
      if quality != None :
         self.quality[iseq]=quality
      return interval
   def col_slice(self,col1,col2) :
      new=self._make_copy_of_all_but_H()
      try :
         new.H = self.H[:,col1:(col2+1)]
         new.base=self.base[col1:(col2+1)]
         new._recompute_base()
         return new
      except :
         return None
   def col_select(self,idx) :
      new=self._make_copy_of_all_but_H()
      try :
         new.H = self.H[:,idx]
         new.base=self.base[idx]
         new._recompute_base()
         return new
      except :
         return None
   def get_histogram_row(self,iseq) :
      return self.H[iseq]
   def __getitem__(self,*arg) :
      if len(arg) == 1 :
         return self.H[arg[0]]
      if len(arg) == 2 :
         return self.H[arg[0],arg[1]]
      return None
   #def __slice__(self,i1,i2) :
      #return self.row_slice(i1,i2)
   def row_slice(self,row1,row2) :
      new=self._make_copy_of_scalars()
      try :
         new.quality = self.quality[row1:(row2+1)]
         new.H = self.H[row1:(row2+1)]
         new.n_seq = len(new.H)
         if self.rows != None :
            new.row=self.rows[row1:(row2+1)]
         else :
            new.rows=None
         return new
      except :
         return None
   def row_select(self,idx) :
      new=self._make_copy_of_scalars()
      try :
         new.quality = self.quality[idx]
         new.H = self.H[idx]
         new.n_seq = len(new.H)
         if self.rows != None :
            new.rows=self.rows[idx]
         else :
            new.rows=None
         if self.rows_delta != None :
            new.rows_delta=self.rows_delta[idx]
         else :
            new.rows_delta=None
         if self.stable_fraction != None :
            new.stable_fraction=self.stable_fraction[idx]
         else :
            new.stable_fraction=None
         return new
      except :
         return None
   def row_number(self,iseq) :
      if self.kind=='differential' :
         return self.H[iseq].sum()
      else :
         return self.H[iseq].max()
   def row_min_max(self,iseq) :
      import numpy as np
      idx=np.where(self.H[iseq]>0)[0]
      if len(idx)==0 :
         return np.array([np.nan,np.nan])
      else :
         imin=idx.min()
         imax=idx.max()
         return np.array([self.base[imin]-self.step/2,self.base[imax]+self.step/2])         
   def finite_base_idx(self) :
      import numpy as np
      return np.where(np.isfinite(self.base))
   def row_difference(self,irow1,irow2,argabsmax=False,absmax=False,exclude_inf=False) :
      import numpy as np
      try :
         diff = np.array(self.H[irow2],dtype='float')-np.array(self.H[irow1],dtype='float')
         if exclude_inf :
            idx = self.finite_base_idx()
            if absmax : 
               return np.abs(diff[idx]).max()
            if argabsmax : 
               return idx[np.abs(diff[idx]).argmax()]
            return diff[idx]
         else :
            if absmax : 
               return np.abs(diff).max()
            if argabsmax : 
               return np.abs(diff).argmax()
            return diff
      except :
         return None 
   def row_sum2(self,irow1,irow2,exclude_inf=False) :
      try :
         if exclude_inf :
            idx = self.finite_base_idx()
            return (self.H[irow1]+self.H[irow2])[idx]
         else :
            return self.H[irow1]+self.H[irow2]
      except :
         return None 
   def row_Neff2(self,irow1,irow2,exclude_inf=False) :
      if not self.Iam['differential'] :
         try :
            return 1./(1./self.Ntot[irow2]+1./self.Ntot[irow1])
         except :
            return None 
         
      try :
         if exclude_inf :
            idx = self.finite_base_idx()
            return 1./((1./(self.H[irow1][idx]).sum()+1./(self.H[irow2])[idx]).sum())
         else :
            return 1./(1./self.H[irow1].sum()+1./self.H[irow2].sum())
      except :
         return None 
      
   def ks_cdf(self,t,mmax=50):#,eps1 = 0.001,eps2 = 1.e-8) :
      """ Kolomogorov Smirnov CDF that, by chance, in the H0 hypothesis that 
          two distributions are equal sqrt(N)*D > t.
          with 
          For one queue distribution
              D = MAX(ABS(Cumulant(X)/N - F(X))
              N = MAX(Cumulant(X)) 
              Cumulant(X) number of samples with value <=X
         For two sample distributions named A and B
              D = MAX(ABS(CumulantA(X)/NA - CumulantA(X)/NA)
              NA = MAX(CumulantA(X)),NB = MAX(CumulantB(X))
              N=1/(1/NA+1/NA) 
         Limit values to accept the H0 hypothesis at ks_cdf< p are
             p     t_max(p)
            ---------------
             0.01   1.63
             0.05   1.36
             
         i.e. if t < t_max(0.01) H0 is valid 
               (the two distributions are equal with 1% probability of false allarm)
         i.e. if t > t_max(0.01) H0 is rejected 
               (the two distributions are not equal with 1% probability of false allarm)
      """
      import math
      #Stop if current term less than EPS1 times previous term
      #Stop if current term changes output by factor less than EPS2
      acc = 0.
      sgn=1.
      for i in range(2*mmax,0,-1) :
         sgn=(-1)**(i-1)
         dacc=sgn*np.exp(-2.*(t*float(i))**2)
         acc+=dacc
      return 2.*acc
   def row_Dvalue2(self,irow1,irow2,exclude_inf=False,bad=None) :
      """returns the 2 queues Dvalue (kolmogorov smirnov D value) for histogram
       rows irow1, irow2
       See ks_cdf for more details.
       
       Returns:
         (D,D*sqrt(N),1.63/sqrt(N),1.36/sqrt(N),N,ks_cdf(D*sqrt(N)))
      """
      import numpy as np
      if not self.Iam['frequency'] and self.Iam['differential'] :
         print "Error: Dvalues are meaningfull just for cumulant and frequency"
         return None
      #################
      #def probks (N_eff,D) :
         #print N_eff,D
         #import math
         #eps1 = 0.001    #Stop if current term less than EPS1 times previous term
         #eps2 = 1.e-8    #Stop if current term changes output by factor less than EPS2
         #en = float(N_eff)**0.5
         #if N_eff == 0 :
            #return 1.
         #Lambda = (en + 0.12 + 0.11/en)*float(D)
         #a2 = -2.*Lambda**2
         #probks = 0.
         #termbf = 0.
         #sign = 1.
         #for j in range(1,101) :
            #term = sign*2*math.exp(a2*float(j)**2)
            #probks = probks + term
            #if ( abs(term) <= eps1*termbf ) or ( abs(term) <= eps2*probks ) :
               #return probks
            #sign = -sign                  #Series alternates in sign
            #termbf = abs(term)
         #probks = np.nan #Sum did not converge after 100 iterations
      #################
      ## va' calcolata sui cumulanti
      D = abs(self.row_difference(irow1,irow2,exclude_inf=False)).max()
      n1=float(self.Ntot[irow1])
      n2=float(self.Ntot[irow2])
      if n1+n2 == 0. : return [D,0.,0.,0.,1.]
      neff = 1/(1/n1+1/n2)
      if D == 0. : return [D,neff,0.,0.,1.]
      sqrtNeff=neff**0.5
      return [D,D*sqrtNeff,1.63/sqrtNeff,1.36/sqrtNeff,neff,self.ks_cdf(D*sqrtNeff)]
   def row_DvalueF(self,irow,F,N):
      #compares the row irow with the cumulant F, numerosity N
      D=abs(self.H[irow]-F).max()
      if N==0:
         neff=float(self.Ntot[irow])
      else:
         n1=float(self.Ntot[irow])
         neff=1/(1/n1+1/float(N))
      if D == 0. : return [D,neff,0.,0.,1.]
      sqrtNeff=neff**0.5
      return [D,D*sqrtNeff,1.63/sqrtNeff,1.36/sqrtNeff,neff,self.ks_cdf(D*sqrtNeff)]
   def row_chisqF(self,irow,F,specific=False):
      #compares the row irow with the cumulant F, numerosity N
      import numpy as np
      D=(np.array(self.H[irow],dtype='float')-F)**2
      N=np.array(self.H[irow],dtype='float')+F
      idx=np.where(N>0)[0]
      if len(idx)==0: return [0.,0.]
      if specific :
         return [(D[idx]/N[idx]).sum()/float(N.sum()),N.sum()]
      return [(D[idx]/N[idx]).sum(),N.sum()]
   
   def row_chisq2(self,irow1,irow2,exclude_inf=False,bad=None) :
      import numpy as np
      diff2 = (self.row_difference(irow1,irow2,exclude_inf=exclude_inf))**2
      sum = self.row_sum2(irow1,irow2,exclude_inf=exclude_inf) 
      idx = np.where(sum > 0)[0]
      if len(idx) == 0 :
         if bad == None :
            return np.nan
         else :
            return np.nan
   def row_mean(self,iseq) :
      num=(self.base[1:-1]*self.H[iseq][1:-1]).sum()
      den=(self.H[iseq][1:-1]).sum()
      return num/float(den) 
   def row_var(self,iseq) :
      dx=(self.base[1:-1]-self.row_mean(iseq))**2
      num=(dx*self.H[iseq][1:-1]).sum()
      den=float((self.H[iseq][1:-1]).sum())
      return num/den
   def row_central_m3(self,iseq) :
    dx=(self.base[1:-1]-self.row_mean(iseq))**3
    num=(dx*self.H[iseq][1:-1]).sum()
    den=float((self.H[iseq][1:-1]).sum())
    return num/den
   def row_central_m4(self,iseq) :
    dx=(self.base[1:-1]-self.row_mean(iseq))**4
    num=(dx*self.H[iseq][1:-1]).sum()
    den=float((self.H[iseq][1:-1]).sum())
    return num/den
   def row_central_mn(self,iseq,n) :
    dx=(self.base[1:-1]-self.row_mean(iseq))**n
    num=(dx*self.H[iseq][1:-1]).sum()
    den=float((self.H[iseq][1:-1]).sum())
    return num/den 
   def row_not_central_mn(self,iseq,n) :
    dx=(self.base[1:-1])**n
    num=(dx*self.H[iseq][1:-1]).sum()
    den=float((self.H[iseq][1:-1]).sum())
    return num/den 
   def row_entropy(self,iseq) :
    import numpy as np
    N=self.H[iseq][1:-1].sum()
    if N==0 :
      return 0.
    freq=self.H[iseq][1:-1]/float(N)
    idx=np.where(self.H[iseq][1:-1]>0)[0]
    return -(freq[idx]*np.log(freq[idx])).sum()/np.log(2.) 
   def min_max(self) :
    import numpy as np
    x=np.zeros([self.n_seq,2])
    for i_row in range(self.n_seq) :
      x[i_row,:]=self.row_min_max(i_row)
    return x.transpose()
   def number(self) :
    x=np.zeros(self.n_seq)
    for i_row in range(self.n_seq) :
      x[i_row]=self.row_number(i_row)
    return x
   def mean(self) :
    import numpy as np
    x=np.zeros(self.n_seq)
    for i_row in range(self.n_seq) :
      x[i_row]=self.row_mean(i_row)
    return x
   def var(self) :
    import numpy as np
    x=np.zeros(self.n_seq)
    for i_row in range(self.n_seq) :
      x[i_row]=self.row_var(i_row)
    return x
   def rms(self) :
    import numpy as np
    return self.var()**0.5
   def central_m3(self) :
    import numpy as np
    x=np.zeros(self.n_seq)
    for i_row in range(self.n_seq) :
      x[i_row]=self.row_central_m3(i_row)
    return x
   def central_m4(self) :
    import numpy as np
    x=np.zeros(self.n_seq)
    for i_row in range(self.n_seq) :
      x[i_row]=self.row_central_m4(i_row)
    return x
   def central_mn(self,n) :
    import numpy as np
    x=np.zeros(self.n_seq)
    for i_row in range(self.n_seq) :
      x[i_row]=self.row_central_mn(i_row,n)
    return x
   def not_central_mn(self,n) :
    import numpy as np
    x=np.zeros(self.n_seq)
    for i_row in range(self.n_seq) :
      x[i_row]=self.row_not_central_mn(i_row,n)
    return x
   def imshow(self,colorbar=True,origin='upper',aspect=None,interpolation='nearest',decimate=None) :
      import numpy as np
      try :
         from matplotlib import pyplot as plt
      except :
         return None
      try :
         dstep=int(decimate)
         idx = np.arange(len(self)/dstep)*dstep
      except :
         dstep=1
         idx = np.arange(len(self))
      plt.imshow(self.H[idx],origin=origin,aspect=aspect,interpolation=interpolation)
      plt.xticks(np.arange(len(self.base)),np.array(self.base,dtype='string'))
      if self.name != None : plt.title(self.name)
      if self.row_name != None : plt.ylabel(self.row_name)
      unit = ''
      if self.col_unit != '' and self.col_unit != None : unit=" [%s] "%self.col_unit
      if self.col_name != None : plt.xlabel(self.col_name+unit)
      plt.plot(0.5*np.ones(2),[-0.5,len(idx)-1+0.5],'w',linewidth=2)
      plt.plot((len(self.base)-2+0.5)*np.ones(2),[-0.5,len(idx)-1+0.5],'w',linewidth=2)
      plt.axis([-0.5,len(self.base)-1+0.5,-0.5,len(idx)-1+0.5])
      if colorbar : plt.colorbar()
   def plot_col(self,iel,iselect=[],xaxis='',symbol='.',label='#') :
      import numpy as np
      try :
         from matplotlib import pyplot as plt
      except :
         return None
      try :
         _label=label
         if _label=='#' : _label = str(iel)
         if _label=='@' : _label = str(self.base[iel])
      except :
         _label=label
      if len(iselect) == 0 :
         try :
            return plt.plot(self.__dict__[xaxis],self.H[:,iel],symbol,label=_label)
         except :
            try :
               return plt.plot(self.H[:,iel],symbol,label=_label)
            except :
               return None
      else :
         try :
            return plt.plot(self.__dict__[xaxis][iselect],self.H[iselect,iel],symbol,label=_label)
         except :
            try :
               return plt.plot(self.H[iselect,iel],symbol,label=_label)
            except :
               return None
   def plot_row(self,iel,iselect=[],xaxis='base',symbol='.',label='#') :
      import numpy as np
      try :
         from matplotlib import pyplot as plt
      except :
         return None
      try :
         _label=label
         if _label=='#' : _label = str(iel)
         if _label=='@' : _label = str(self.rows[iel])
      except :
         _label=label
      if len(iselect) == 0 :
         try :
            return plt.plot(self.__dict__[xaxis],self.H[iel],symbol,label=_label)
         except :
            try :
               return plt.plot(self.H[iel],symbol,label=_label)
            except :
               return None
      else :
         try :
            return plt.plot(self.__dict__[xaxis][iselect],self.H[iselect,iel],symbol,label=_label)
         except :
            try :
               return plt.plot(self.H[iselect,iel],symbol,_label)
            except :
               return None
   def entropy(self) :
    import numpy as np
    x=np.zeros(self.n_seq)
    for i_row in range(self.n_seq) :
      x[i_row]=self.row_entropy(i_row)
    return x
   def gaussian_entropy(self,sheppardCorrection=False) :
    import numpy as np
    sigmaq = self.var() 
    if sheppardCorrection :
      sigmaq -= self.step**2/12.
    sigmaq=sigmaq**0.5/self.step
    return np.log((2.*np.pi*np.exp(1.))**0.5*sigmaq)/np.log(2.)
   def cumulant(self) :
      "extract the cumulant neglecting the lower and upper limit of base"
      import copy
      new=copy.deepcopy(self)
      if new.Iam['differential'] :
         if new.Iam['frequency'] :
            new.kind='cumulant:frequency'
         else :
            new.kind='cumulant:number'
         new.Iam['differential']=False
         Ht=new.H.transpose()
         if not self.Iam['frequency'] and abs(new.base[0])==np.inf : new.Ninf = copy.deepcopy(Ht[0])
         if not self.Iam['frequency'] and abs(new.base[-1])==np.inf : new.Nsup = copy.deepcopy(Ht[-1])
         inf = 1 if abs(new.base[0])==np.inf else 0 
         sup = len(new.base) if abs(new.base[-1])==np.inf else len(new.base)+1
         for i in range(inf,sup) : Ht[i]+=Ht[i-1]
         if not self.Iam['frequency'] : new.Ntot=copy.deepcopy(Ht[-1])
         new.H = Ht.transpose()
      return new
   def frequency(self) :
      "convert to frequency neglecting lower and upper limit of base"
      import copy
      new=copy.deepcopy(self)
      if not new.Iam['frequency'] :
         if new.Iam['differential'] :
            new.kind='differential:frequency'
         else :
            new.kind='cumulant:frequency'
         new.Iam['frequency']=True
         new.H = np.array(new.H,dtype='float')
         new.celltype='float'
         if self.Iam['differential'] and abs(new.base[0])==np.inf and new.Ninf == None : new.Ninf = copy.deepcopy(new.H[:,0])
         if self.Iam['differential'] and abs(new.base[-1])==np.inf and new.Nsup == None : new.Nsup = copy.deepcopy(new.H[:,-1])
         inf = 1 if abs(new.base[0])==np.inf else 0 
         sup = len(new.base) if abs(new.base[-1])==np.inf else len(new.base)+1
         if self.Iam['differential'] and new.Ntot == None : 
            new.Ntot=np.zeros(len(new))
            for i in range(0,len(new)) : 
               new.Ntot[i] = float(self.H[i,inf:sup].sum())
         for i in range(0,len(new)) : new.H[i]*=1./float(new.Ntot[i])
      return new
   def coadd_rows(self,step) :
    "return an object with coadded rows to reduce the number of rows"
    import numpy as np
    if step <= 0 : return None
    n_seq = int(np.ceil(float(self.n_seq)/float(step)))
    new = HistogramSequence(n_seq,self.base,n_over=self.n_over,celltype=self.celltype,name=self.name)
    if self.rows != None :
       new.rows=np.zeros(new.n_seq)
    if self.heights != None :
       new.rows=np.zeros(new.n_seq)
    if self.stable_fraction != None :
       new.stable_fraction=np.zeros(new.n_seq)
    for i_new_seq in range(new.n_seq) :
       if (i_new_seq+1)*step-1 < len(self) :
          iseq=np.arange(i_new_seq*step,(i_new_seq+1)*step)
       else :
          iseq=np.arange(i_new_seq*step,len(self))
       if self.rows != None :
          rmin = self.rows[iseq].min()
          rmax = self.rows[iseq].max()
          new.rows[i_new_seq]=(rmax+rmin)/2.
       if self.heights != None :
          new.heights[i_new_seq]=self.heights[iseq].sum()
       for i_iseq in iseq :
          new.H[i_new_seq]+=self.H[i_iseq]
       if self.stable_fraction != None :
          new.stable_fraction[i_new_seq]=(self.stable_fraction[iseq]*self.rows_delta[iseq]).sum()/float(self.rows_delta[iseq].sum())
    return new
   def get_rows(self,*arg) :
    "return an object with a subset of rows"
    if len(arg) == 0 or len(arg) > 2: return None
    row1=arg[0]
    try : 
      row2=arg[1]
    except :
      row2=arg[0]
    if row1 < 0 or row2 < 0 : return None
    if row1 >= len(self) or row2 >= len(self) : return None
    if row1 > row2 : return None
    n_seq = row2-row1+1
    new = HistogramSequence(n_seq,self.base,n_over=self.n_over)
    new.H=self.H[row1:(row2+1)]
    if self.heights != None : new.heights=self.heights[row1:(row2+1)]
    if self.rows != None : new.rows=self.rows[row1:(row2+1)]
    return new
   def append(self,Hin) :
      "append an input histogram"
      import numpy as np
      for k in self._keys_vectors() :
         print k,
         if self.__dict__.has_key(k) and Hin.__dict__.has_key(k) :
            print 'haS',
            if self.__dict__[k]!=None and Hin.__dict__[k]!=None :
               print 'good'
               self.__dict__[k] = np.concatenate((self.__dict__[k],Hin.__dict__[k]),axis=0)
            else :
               print 'none',
         else :
            print 'not has'
      self.n_seq+=Hin.n_seq

   
import numpy as np
if __name__=='__main__' :
  x=np.array([96.99,105.01])
  HS=HistogramSequence(5,100.,102.,1)
  print HS.fill(0,np.array([96,99,101,100]))
  print HS.fill(1,np.array([98,99,100,101,102,103,104]))
  print HS.fill(2,np.array([95.99,99.9,100.02]))
  print HS.fill(3,np.array([99.98,102.11,102.21]))
  print HS.fill(4,np.array([101.97,104.12]))
  print HS.H
  
  nseq=8
  x=np.array([])
  HS=HistogramSequence(nseq,0.,3.,1,col_name='X',col_unit='arbitrary',row_name='row')
  for i in range(nseq) :
    HS.fill(i,np.repeat([-2,-1,0.,1,2,3,4,5],i+1))
