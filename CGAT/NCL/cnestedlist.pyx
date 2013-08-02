#cython: embedsignature=True
cimport cython

###############################
# Could not make .pxd file to be found in gpipe/setup.py, so including it here:
# Fields to extension classed have been added to each class.
###############################
cdef extern from "string.h":
  ctypedef int size_t
  void *memcpy(void *dst,void *src,size_t len)
  void *memmove(void *dst,void *src,size_t len)
  void *memset(void *b,int c,size_t len)

cdef extern from "stdlib.h":
  void free(void *)
  void *malloc(size_t)
  void *calloc(size_t,size_t)
  void *realloc(void *,size_t)
  int c_abs "abs" (int)
  void qsort(void *base, size_t nmemb, size_t size,
             int (*compar)(void *,void *))

cdef extern from "stdio.h":
  ctypedef struct FILE:
    pass
  FILE *fopen(char *,char *)
  int fclose(FILE *)
  int sscanf(char *str,char *fmt,...)
  int sprintf(char *str,char *fmt,...)
  int fprintf(FILE *ifile,char *fmt,...)
  char *fgets(char *str,int size,FILE *ifile)

cdef extern from "string.h":
  int strcmp(char *s1, char *s2)
  int strncmp(char *s1,char *s2,size_t len)
  char *strcpy(char *dest,char *src)
  char *strdup(char *)
  char *strcat(char *,char *)

cdef extern from "intervaldb.h":
  ctypedef struct IntervalMap:
    int start
    int end
    int target_id
    int sublist

  ctypedef struct IntervalIndex:
    int start
    int end

  ctypedef struct SublistHeader:
    int start
    int len

  ctypedef struct SubheaderFile:
    pass
  
  ctypedef struct IntervalDBFile:
    int n
    int ntop
    int nlists
    int div
    int nii
    IntervalIndex *ii
    SublistHeader *subheader
    SubheaderFile subheader_file
    FILE *ifile_idb

  ctypedef struct IntervalIterator:
    pass

  ctypedef struct FilePtrRecord:
    FILE *ifile
    int left
    int right
    int ihead
    char *filename

  int imstart_qsort_cmp(void *void_a,void *void_b)
  int target_qsort_cmp(void *void_a,void *void_b)
  IntervalMap *read_intervals(int n,FILE *ifile) except NULL
  SublistHeader *build_nested_list(IntervalMap im[],int n,int *p_n,int *p_nlists) except NULL
  SublistHeader *build_nested_list_inplace(IntervalMap im[],int n,int *p_n,int *p_nlists) except NULL
  IntervalMap *interval_map_alloc(int n) except NULL
  IntervalIterator *interval_iterator_alloc() except NULL
  int free_interval_iterator(IntervalIterator *it)
  IntervalIterator *reset_interval_iterator(IntervalIterator *it)
  int find_intervals(IntervalIterator *it0,int start,int end,IntervalMap im[],int n,SublistHeader subheader[],int nlists,IntervalMap buf[],int nbuf,int *p_nreturn,IntervalIterator **it_return) except -1
  char *write_binary_files(IntervalMap im[],int n,int ntop,int div,SublistHeader *subheader,int nlists,char filestem[])
  IntervalDBFile *read_binary_files(char filestem[],char err_msg[],int subheader_nblock) except NULL
  int free_interval_dbfile(IntervalDBFile *db_file)
  int find_file_intervals(IntervalIterator *it0,int start,int end,IntervalIndex ii[],int nii,SublistHeader subheader[],int nlists,SubheaderFile *subheader_file,int ntop,int div,FILE *ifile,IntervalMap buf[],int nbuf,int *p_nreturn,IntervalIterator **it_return) except -1
  int write_padded_binary(IntervalMap im[],int n,int div,FILE *ifile)
  int read_imdiv(FILE *ifile,IntervalMap imdiv[],int div,int i_div,int ntop)
  int save_text_file(char filestem[],char basestem[],char err_msg[],FILE *ofile)
  int text_file_to_binaries(FILE *infile,char err_msg[])
  int C_int_max

###############################
### .pxd end
###############################

cdef class IntervalDBIterator:
  """Iterator over intervals from an NCL."""

  cdef IntervalIterator *it,*it_alloc
  cdef IntervalMap im_buf[1024]
  cdef int ihit,nhit,start,end
  cdef IntervalDB db

  def __cinit__(self,int start,int end,IntervalDB db not None):
    self.it=interval_iterator_alloc()
    self.it_alloc=self.it
    self.start=start
    self.end=end
    self.db=db

    # stop if interval is empty
    if self.start >= self.end: 
      raise IndexError( "invalid interval (%i,%i)" % (self.start, self.end) )

  def __iter__(self):
    return self 

  cdef int cnext(self): # C VERSION OF ITERATOR next METHOD RETURNS INDEX

    cdef int i
    if self.ihit>=self.nhit: # TRY TO GET ONE MORE BUFFER CHUNK OF HITS
      if self.it==NULL: # ITERATOR IS EXHAUSTED
        return -1
      find_intervals(self.it,self.start,self.end,self.db.im,self.db.ntop,
                     self.db.subheader,self.db.nlists,self.im_buf,1024,
                     &(self.nhit),&(self.it)) # GET NEXT BUFFER CHUNK
      self.ihit=0 # START ITERATING FROM START OF BUFFER
    if self.ihit<self.nhit: # RETURN NEXT ITEM FROM BUFFER
      i=self.ihit
      self.ihit = self.ihit+1 # ADVANCE THE BUFFER COUNTER
      return i
    else: # BUFFER WAS EMPTY, NO HITS TO ITERATE OVER...
      return -1

  # PYTHON VERSION OF next RETURNS HIT AS A TUPLE
  def __next__(self): # PYREX USES THIS NON-STANDARD NAME INSTEAD OF next()!!!
    cdef int i
    i=self.cnext()
    if i>=0:
      return (self.im_buf[i].start,self.im_buf[i].end,self.im_buf[i].target_id )
    else:
      raise StopIteration

  def __dealloc__(self):
    'remember: dealloc cannot call other methods!'
    free_interval_iterator(self.it_alloc)


cdef class IntervalDB:
  """a NCL database in memory.

  The initial database is empty. Call any of the ``from``
  methods to populate the database.
  """
  cdef int n
  cdef int ntop
  cdef int nlists
  cdef IntervalMap *im
  cdef SublistHeader *subheader


  def __cinit__(self,**kwargs):
    self.n = 0

  def fromsortedfile( self, filename, nsize, **kwargs):
    """build database from *nsize* tuples in *filename*.

    Each line in the file should contain three whitespace separated integers
    with ``start``, ``end`` and ``index`` of an interval.

    The intervals have to be sorted.
    """
    cdef FILE *ifile
    cdef int i
    assert nsize > 0, "number of tuples must be larger than 0"
    self.n=nsize

    ifile=fopen(filename,"r")
    if ifile:
      self.im=read_intervals(self.n,ifile)
      fclose(ifile)
      if self.im!=NULL:
        self.runBuildMethod(**kwargs)
    else:
      raise IOError('could not open file %s' % filename)

  def fromlist(self, l, **kwargs):
    '''build from data in l.

    The data should be given as list of tuple of tuples
    ``(start,end,id)``.

    see :meth:runBuildMethod for *kwargs*.
    '''
    cdef int i
    self.close() # DUMP OUR EXISTING MEMORY 
    self.n=len(l)
    self.im=interval_map_alloc(self.n)
    if self.im==NULL:
      raise MemoryError('unable to allocate IntervalMap[%d]' % self.n)
    i=0
    for t in l:
      self.im[i].start=t[0]
      self.im[i].end=t[1]
      self.im[i].target_id=t[2]
      self.im[i].sublist= -1
      i=i+1
    self.runBuildMethod(**kwargs)

  def runBuildMethod(self, buildInPlace=True):
    '''build either in-place if *buildInPlace == True* or using older build method

    TODO: this method currently does not check for invalid intervals.
    '''
    if buildInPlace:
      self.subheader=build_nested_list_inplace(self.im,self.n,&(self.ntop),&(self.nlists))
    else:
      self.subheader=build_nested_list(self.im,self.n,&(self.ntop),&(self.nlists))

  def fromunsortedfile(self, filename,int n,**kwargs):
    """build database from *nsize* tuples in *filename*.

    The intervals need not be sorted.
    """
    cdef FILE *ifile
    cdef int i
    cdef IntervalMap *im_new
    self.close()
    ifile=fopen(filename,'r')
    if ifile==NULL:
      raise IOError('unable to open '+filename)
    im_new=interval_map_alloc(n)
    if im_new==NULL:
      raise MemoryError('unable to allocate IntervalMap[%d]' % n)
    i=read_imdiv(ifile,im_new,n,0,n)
    fclose(ifile)
    if i!=n:
      raise IOError('IntervalMap file corrupted?')
    self.n=n
    self.im=im_new
    self.runBuildMethod(**kwargs)

  def find_overlap(self,int start,int end):
    """find intervals in database overlapping with *start* and *end*.
    
    returns an :class:`ncl.IntervalDBIterator`
    """
    self.check_nonempty() # RAISE EXCEPTION IF NO DATA
    return IntervalDBIterator(start,end,self)

  def find_overlap_list(self,int start,int end):
    """return list of intervals in database overlapping with *start* and *end*.
    """
    cdef int i,nhit
    cdef IntervalIterator *it,*it_alloc
    cdef IntervalMap im_buf[1024]
    self.check_nonempty() # RAISE EXCEPTION IF NO DATA
    it=interval_iterator_alloc()
    it_alloc=it
    l=[] # LIST OF RESULTS TO HAND BACK
    while it:
      find_intervals(it,start,end,self.im,self.ntop,
                     self.subheader,self.nlists,im_buf,1024,
                     &(nhit),&(it)) # GET NEXT BUFFER CHUNK
      for i from 0 <= i < nhit:
        l.append((im_buf[i].start,im_buf[i].end,im_buf[i].target_id))
    free_interval_iterator(it_alloc)
    return l
        
  def check_nonempty(self):
    """return True if the database is empty."""
    if self.im:
      return True
    else:
      msg='empty IntervalDB, not searchable!'
      raise IndexError(msg)

  def write_binaries(self,filestem,div=256):
    """commit database to filesystem using *filestem* as 
    root name of files.

    *div* refers to the block size.
    """
    cdef char *err_msg
    err_msg=write_binary_files(self.im,self.n,self.ntop,div,
                               self.subheader,self.nlists,filestem)
    if err_msg:
      raise IOError(err_msg)

  def __dealloc__(self):
    'remember: dealloc cannot call other methods!'
    if self.subheader:
      free(self.subheader)
    if self.im:
      free(self.im)
    
  def close(self):
    if self.subheader:
      free(self.subheader)
    if self.im:

      free(self.im)
    self.subheader=NULL
    self.im=NULL

    return None

cdef class IntervalFileDB:
  """a NCL database on the filesystem identified by *filestem*.

  NCL disk-based databases can only be opened for reading.
  create the database by building it in memory using :class:`ncl.IntervalDB` and then saving
  it to disk.
  """
  cdef IntervalDBFile *db

  def __cinit__(self,filestem=None):
    if filestem is not None:
      self.open(filestem)

  def open(self,filestem):
    cdef char err_msg[1024]
    self.db=read_binary_files(filestem,err_msg,1024)
    if self.db==NULL:
      raise IOError(err_msg)

  def find_overlap(self,int start,int end):
    """find intervals in database overlapping with *start* and *end*.
    
    returns an :class:`ncl.IntervalDBIterator`
    """
    self.check_nonempty() # RAISE EXCEPTION IF NO DATA
    return IntervalFileDBIterator(start,end,self)

  def find_overlap_list(self,int start,int end):
    """return list of intervals in database overlapping with *start* and *end*.
    """
    cdef int i,nhit
    cdef IntervalIterator *it,*it_alloc
    cdef IntervalMap im_buf[1024]
    self.check_nonempty() # RAISE EXCEPTION IF NO DATA
    it=interval_iterator_alloc()
    it_alloc=it
    l=[] # LIST OF RESULTS TO HAND BACK
    while it:
      find_file_intervals(it,start,end,self.db[0].ii,self.db[0].nii,
                          self.db[0].subheader,self.db[0].nlists,
                          &(self.db[0].subheader_file),
                          self.db[0].ntop,self.db[0].div,
                          self.db[0].ifile_idb,im_buf,1024,
                          &(nhit),&(it)) # GET NEXT BUFFER CHUNK
      for i from 0 <= i < nhit:
        l.append((im_buf[i].start,im_buf[i].end,im_buf[i].target_id))
    free_interval_iterator(it_alloc)
    return l

  def check_nonempty(self):
    if self.db==NULL:
      raise IndexError('empty IntervalFileDB, not searchable!')

  def close(self):
    if self.db:
      free_interval_dbfile(self.db)
    self.db=NULL

  def __dealloc__(self):
    'remember: dealloc cannot call other methods!'
    if self.db:
      free_interval_dbfile(self.db)


cdef class IntervalFileDBIterator:
  """disk based intervalDB."""

  cdef IntervalIterator *it,*it_alloc
  cdef IntervalMap *im_buf
  cdef int ihit,nhit,start,end,nbuf
  cdef IntervalFileDB db
  cdef IntervalDB idb

  cdef int restart(self,int start,int end,IntervalFileDB db) except -2
  cdef int reset(self) except -2
  cdef int cnext(self,int *pkeep)
  cdef int extend(self,int ikeep)
  cdef int saveInterval(self,int start,int end,int target_id)
  cdef int nextBlock(self,int *pkeep) except -2
  cdef IntervalMap *getIntervalMap(self)
  cdef int loadAll(self) except -1
  cdef int copy(self,IntervalFileDBIterator src)

  def __cinit__(self,
              int start,int end,
              IntervalFileDB db=None,
              int nbuffer=1024,rawIvals=None):
    cdef int i

    if start >= end: 
      raise IndexError( "invalid interval (%i,%i)" % (start, end) )

    ns = None
    self.it=interval_iterator_alloc()
    self.it_alloc=self.it
    self.start=start
    self.end=end
    self.db=db
    if ns is not None:
      if ns.idb is not None:
        self.idb=ns.idb
      elif ns.db is None:
        ns.forceLoad()
      self.db=ns.db
    if rawIvals is not None and len(rawIvals)>nbuffer:
      nbuffer=len(rawIvals)
    self.im_buf=interval_map_alloc(nbuffer)
    self.nbuf=nbuffer
    if rawIvals is not None:
      i=0
      for ival in rawIvals:
        self.im_buf[i].start=ival[0] # SAVE INTERVAL INFO
        self.im_buf[i].end=ival[1]
        self.im_buf[i].target_id=ival[2]
        i=i+1
      self.nhit=i # TOTAL NUMBER OF INTERVALS STORED

  cdef int restart(self,int start,int end,IntervalFileDB db) except -2:
    'reuse this iterator for another search without reallocing memory'
    self.nhit=0 # FLUSH ANY EXISTING DATA
    self.start=start
    self.end=end
    self.db=db
    self.it=self.it_alloc # REUSE OUR CURRENT ITERATOR
    reset_interval_iterator(self.it) # RESET IT FOR REUSE
    return 0

  cdef int reset(self) except -2:
    'flush the buffer so we can reuse this iterator'
    self.nhit=0
    return 0

  cdef int extend(self,int ikeep):
    'expand the buffer if necessary, keeping elements [ikeep:nbuf]'
    cdef int length,istart
    cdef IntervalMap *new_buf
    istart=self.nbuf-ikeep
    length=sizeof(IntervalMap)*istart # #BYTES WE MUST KEEP
    if ikeep==0: # KEEPING THE WHOLE BUFFER, SO MUST ALLOCATE NEW SPACE
      new_buf=<IntervalMap *>realloc(self.im_buf,2*length) # DOUBLE OUR BUFFER
      if new_buf==NULL:
        raise MemoryError('out of memory')
      self.im_buf=new_buf
      self.nbuf=2*self.nbuf
    elif ikeep<8: # RUNNING OUT OF ROOM, SO EXPAND BUFFER
      self.nbuf=2*self.nbuf
      new_buf=interval_map_alloc(self.nbuf)
      memcpy(new_buf,self.im_buf+ikeep,length)
      free(self.im_buf)
      self.im_buf=new_buf
    else: # JUST SHIFT [ikeep:] SLICE OF BUFFER TO FRONT [0:istart]
      memmove(self.im_buf,self.im_buf+ikeep,length)
    return istart # RETURN START OF EMPTY BLOCK WHERE WE CAN ADD NEW DATA

  cdef int saveInterval(self,int start,int end,int target_id):
    'save an interval, expanding array if necessary'
    cdef int i
    if self.nhit>=self.nbuf: # EXPAND ARRAY IF NECESSARY
      self.extend(0)
    i=self.nhit
    self.im_buf[i].start=start # SAVE INTERVAL INFO
    self.im_buf[i].end=end
    self.im_buf[i].target_id=target_id
    self.nhit = i+1
    return self.nhit

  cdef int nextBlock(self,int *pkeep) except -2:
    'load one more block of overlapping intervals'
    cdef int i
    if self.it==NULL: # ITERATOR IS EXHAUSTED
      return -1
    if pkeep and pkeep[0]>=0 and pkeep[0]<self.nhit: #MUST KEEP [ikeep:] SLICE
      i=self.extend(pkeep[0]) # MOVE SLICE TO THE FRONT
    else: # WE CAN USE THE WHOLE BUFFER
      i=0
    if self.db is not None: # ON-DISK DATABASE
      find_file_intervals(self.it,self.start,self.end,
                          self.db.db[0].ii,self.db.db[0].nii,
                          self.db.db[0].subheader,self.db.db[0].nlists,
                          &(self.db.db[0].subheader_file),
                          self.db.db[0].ntop,self.db.db[0].div,
                          self.db.db[0].ifile_idb,
                          self.im_buf+i,self.nbuf-i,
                          &(self.nhit),&(self.it)) # GET NEXT BUFFER CHUNK
    elif self.idb is not None: # IN-MEMORY DATABASE
      find_intervals(self.it,self.start,self.end,self.idb.im,self.idb.ntop,
                     self.idb.subheader,self.idb.nlists,self.im_buf+i,self.nbuf-i,
                     &(self.nhit),&(self.it)) # GET NEXT BUFFER CHUNK
    else:
      raise IOError('Iterator has no database!  Please provide a db argument.')
    self.nhit=self.nhit+i # TOTAL #HITS IN THE BUFFER
    self.ihit=i # START ITERATING FROM START OF NEW HITS
    if pkeep and pkeep[0]>=0: # RESET ikeep INDEX TO START OF BUFFER
      pkeep[0]=0
    return self.nhit-self.ihit # RETURN #NEW HITS IN NEXT BLOCK

  cdef IntervalMap *getIntervalMap(self):
    '''return the IntervalMap array loaded by iterator,
    and release it from iterator.  User must free the array!'''
    cdef int len
    cdef IntervalMap *im
    if self.nhit==0: # NO HITS
      return NULL
    elif self.nhit<self.nbuf: # LARGER BUFFER THAN WE ACTUALLY NEED
      len=sizeof(IntervalMap)*self.nhit # COMPUTE FINAL SIZE
      im=<IntervalMap *>realloc(self.im_buf,len) # COMPACT TO FINAL SIZE
    else: # JUST HAND BACK OUR FULL BUFFER
      im=self.im_buf
    self.im_buf=NULL # RELEASE THIS STORAGE FROM ITERATOR; USER MUST FREE IT!
    return im # HAND BACK THE STORAGE

  cdef int loadAll(self) except -1:
    'load all overlapping interval hits, return count of hits'
    cdef int len,ikeep
    len=1
    ikeep=0 # DON'T LET extend DISCARD ANY HITS, KEEP THEM ALL!
    while len>0: # LOAD BLOCKS UNTIL NO MORE...
      len=self.nextBlock(&ikeep) # LOAD ANOTHER BLOCK OF INTERVALS
    return self.nhit

  cdef int cnext(self,int *pkeep): # C VERSION OF ITERATOR next METHOD
    'get one more overlapping interval'

    cdef int i
    if self.ihit>=self.nhit: # TRY TO GET ONE MORE BUFFER CHUNK OF HITS
      self.nextBlock(pkeep) # LOAD THE NEXT BLOCK IF ANY
    if self.ihit<self.nhit: # RETURN NEXT ITEM FROM BUFFER
      i=self.ihit
      self.ihit = self.ihit+1 # ADVANCE THE BUFFER COUNTER
      return i
    else: # BUFFER WAS EMPTY, NO HITS TO ITERATE OVER...
      return -1

  cdef int copy(self,IntervalFileDBIterator src):
    'copy items from src to this iterator buffer'
    cdef IntervalMap *new_buf
    if src is None:
      raise ValueError('src is None!  Debug!!')
    if src.nhit>self.nbuf: # NEED TO EXPAND OUR BUFFER
      new_buf=<IntervalMap *>realloc(self.im_buf,src.nhit*sizeof(IntervalMap))
      if new_buf==NULL:
        raise MemoryError('out of memory')
      self.im_buf=new_buf # RECORD NEW BUFFER LOCATION AND SIZE
      self.nbuf=src.nhit
    self.nhit=src.nhit # COPY ARRAY AND SET CORRECT SIZE
    if src.nhit>0: # ONLY COPY IF NON-EMPTY
      memcpy(self.im_buf,src.im_buf,src.nhit*sizeof(IntervalMap))
    return 0

  def __iter__(self):
    return self

  # PYTHON VERSION OF next RETURNS HIT AS A TUPLE
  def __next__(self): # PYREX USES THIS NON-STANDARD NAME INSTEAD OF next()!!!
    cdef int i
    i=self.cnext(NULL)
    if i>=0:
      return (self.im_buf[i].start,self.im_buf[i].end,self.im_buf[i].target_id )
    else:
      raise StopIteration

  def __dealloc__(self):
    'remember: dealloc cannot call other methods!'
    free_interval_iterator(self.it_alloc)
    if self.im_buf:
      free(self.im_buf)

      



