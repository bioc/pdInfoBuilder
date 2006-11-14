#!/usr/bin/env python
"""
5 April, 2006, Seth:
    Rename to readAffySnpCsv.py and reformat long lines.

March 18, ross: 
    converted to bioc standardish columns for Seth's affyanno R
    code also added mysql option since that's what I'm using for
    my own (non-R) projects eeek
    http://www.sqlite.org/lang_createindex.html says sqlite
    recreates indexes on opening a database - if this takes time,
    it might be very bad.. No, it's not a problem. Each index adds
    10MB or so but opens fast Added BadHeaderError if any problem
    finding needed column names - tested lightly

ross lazarus me fecit march 13 2006:

* read affy 500k anno files into sqlite

* psyco makes a small (about 5% faster) difference if you can load
  it

* trots along at around 3500 rows/second for loading on non-nfs
  dual athlon 1800GHZ starts faster (5300 rows/sec) on my windoze
  box but drops to around the same load speed as the file grows

* genotypes might be better read into one table/chromosome since
  sqlite slows down badly with large files - much more than mysql
  eg.  but we probably need one mongo table (tm) for annotation
  since we won't always know the chromosome ahead of time although
  we could build a mapping table just for that?

* the associated gene field and the microsatellite fields contain
  an arbitrary number of references separated by forward slash
  pairs and triplets probably need to make separate reference
  tables for these?  first pass does not do that yet..

TODO mar13 2006:

   add indices - but if autoreading column headers, need to think
   about how to do it safely if affy changes those names

   add separate tables for genes and microsatellites - they are
   multiple entries per row with fields separated by "//"

   and entries separated by "///"
"""
import sys, glob, time
import logging

logging.basicConfig(level=logging.DEBUG)
LOG = logging.getLogger("readAffyAnno")


try:
    import psyco
    # adds about 5-10% performance on the machines I tested
    psyco.full() 
except:
    LOG.info('psyco not available')

debug = True


class completeaffyanno:
    """simple class to handle affy 500k mapping annotation files as sqlite
    tables This version reads ALL headers into a sqlite table There's
    a much uglier version which produces agreed column names including
    some derived ones, and optionally uses MySQL instead of sqlite
    below

    """

    class affyannoF:
        """iterable for csv elements to return csv file lines as field lists
        sqlite expects for executemany Inserts a Null value for
        primary key which sqlite will autoincrement

        """
        def __init__(self,fname=''):
            self.f = file(fname,'r')
            self.header = self.f.next()
            self.lnum = 1
            self.fname = fname
            self.started = time.time()

        def next(self):
            try:
                line = self.f.next()
                self.lnum += 1
            except:
                raise StopIteration
            ll = line.strip().split('","')
            ll = [x.replace('"','') for x in ll] # remove quotes from csv
            ll.insert(0,None) # for autoincrement pk
            if self.lnum % 10000 == 0:
                if debug:
                    dur = time.time() - self.started
                    print ('affyannoF reading line %d of %s = %5.2f recs/sec' 
                           % (self.lnum,self.fname,self.lnum/dur))
            return ll

        def __iter__(self):
            return self
            
    def __init__(self,dbname='affy',tablename='SNP500k',annoflist = []):
        """
        """
        self.dbname = dbname
        self.tablename = tablename
        self.annoflist = annoflist
        # open annotation file as an iterable that returns a decsved
        # field list
        f = affyannoF(annoflist[0]) 
        affyannohead = f.header # freshly read
        # is csv
        fieldnames = affyannohead.strip().split(',') 
        # spaces won't work in column names
        fieldnames = [x.replace(' ','_') for x in fieldnames]
        # these are illegal 
        fieldnames = [x.replace('-','') for x in fieldnames] 
        # get rid of quotes
        fieldnames = [x.replace('"','') for x in fieldnames] 
        # remember to insert a null here always
        fieldnames.insert(0,'ID') 
        nfields = len(fieldnames)
        # leave room for pk
        fieldtypes = ['text']*(nfields-1) 
        # make an autoincrement if
        fieldtypes.insert(0,'INTEGER PRIMARY KEY NOT NULL') 
        self.fieldnames = fieldnames
        self.fieldtypes = fieldtypes
        self.createvars = ','.join(['%s %s' % (self.fieldnames[x],self.fieldtypes[x]) for x in range(nfields)])
        # sqlite uses qmarks
        varplaces = ','.join(['?']*nfields) 
        self.insertsql = ('INSERT INTO %s (%s) values (%s)' 
                          % (self.tablename,','.join(self.fieldnames),varplaces))
        if debug:
            LOG.debug('AffyAnnoF initialized')
            LOG.debug('Got fieldnames = %s' % fieldnames)
            LOG.debug('createvars = %s' % self.createvars)
            LOG.debug('insertsql = %s' % self.insertsql)
            
    def build(self):
        """create and fill a sqlite db from the mapping annotation files
        downloaded from
        http://www.affymetrix.com/support/technical/byproduct.affx?product=500k

        """
        # Create a connection to the database file 
        con = sqlite.connect(self.dbname)
        # Get a Cursor object that operates in the context of Connection con:
        cur = con.cursor()
        try:
            sql = 'drop table %s' % self.tablename
            cur.execute(sql)
        except:
            if debug:
                LOG.debug('no %s to drop this time!' % self.tablename)
            else:
                pass
        sql = 'create table %s (%s);' % (self.tablename,self.createvars)
        cur.execute(sql)
        for fname in self.annoflist:
            f = affyannoF(fname)
            # just pass a file iterator!@
            cur.executemany(self.insertsql, f) 
            con.commit() # only at end of each file to speed things up
        cur.close()
        con.close()


class BadHeaderError(Exception):
    """catch bad headers
    """
    pass

      
class affyanno:
    """class to handle affy 500k mapping annotation files and write as
    sqlite tables This version reads agreed column names including
    some derived ones, and optionally uses MySQL instead of sqlite
    below 

    """
 
    class affyannoF:
        """iterable for csv elements to return csv file lines as field lists
        mysql and sqlite expect for executemany Fugly hack of what
        was elegant necessitated by need to translate some column
        names and construct some derived ones...guido, forgive me
        for I have sinned..  Optionally inserts a Null value for
        primary key which sqlite will autoincrement affy csv file
        header: "Probe Set ID","Affy SNP ID","dbSNP RS
        ID","Chromosome","Genome Version", "DB SNP
        Version","Physical Position","Strand","ChrX
        pseudo-autosomal region", "Cytoband","Flank","Allele
        A","Allele B","Associated Gene","Genetic
        Map","Microsatellite", "Fragment Length Start Stop","Freq
        Asian","Freq AfAm","Freq Cauc","Het Asian","Het AfAm",
       "Het Cauc","Num chrm Asian","Num chrm AfAm","Num chrm Cauc"
        NOTE: this version ONLY returns probe_set_id, affy_SNP_id, dbsnpid, 
        chrom, pos, strand as defined in wewant below. Column names from 
        current version of
        http://wiki.fhcrc.org/bioc/SNP_Annotation_Package_Discussion
        probeSetId, affySnpId, dbSnpId, chrom, start, end, alleles, strand, 
        flank

        """

        def __init__(self,fname='', addnullpk=0):
            # this determines which fields we keep!
            self.wewant = ['Probe_Set_ID', 'Affy_SNP_ID', 'dbSNP_RS_ID', 
                           'Chromosome', 'Physical_Position','Strand',
                           'Allele_A', 'Allele_B']
            # standard column names: for now just lowercase
            self.colnames = [ x.lower() for x in self.wewant ]

            self.f = file(fname,'r')
            header = self.f.next()
            # is csv
            fieldnames = header.strip().split(',') 
            # cannot have spaces..
            fieldnames = [x.replace(' ','_') for x in fieldnames]
            # these are not good in sql column names
            fieldnames = [x.replace('-','') for x in fieldnames] 
            # get rid of quotes
            fieldnames = [x.replace('"','') for x in fieldnames] 
            if debug:
                LOG.debug('header=%s' % fieldnames)
                LOG.debug('wewant=%s' % self.wewant)
            try: # make indexes for fields wewant
                # keep these fields ONLY
                self.iwewant = [fieldnames.index(x) for x in self.wewant] 
            except:
                raise BadHeaderError(fname)
            self.lnum = 1
            self.fname = fname
            self.started = time.time()
            # flag for autoincrement pk
            self.addnullpk = addnullpk 

        def next(self):
            try:
                line = self.f.next()
                self.lnum += 1
            except:
                raise StopIteration
            ll = line.strip().split('","')
            # remove remaining quotes from csv
            ll = [x.replace('"','') for x in ll]
            # translate missing value indicator
            def missingToNone(x):
                if x != '---':
                    return x
                return None
            ll = [ missingToNone(x) for x in ll ]
            # filter only the fields we want to save
            affyrec = [ll[x] for x in self.iwewant]
            if self.addnullpk:
                ll.insert(0, None) # for autoincrement pk
            if self.lnum % 10000 == 0:
                if debug:
                    dur = time.time() - self.started
                    print ('affyannoF reading line %d of %s = %5.2f recs/sec' 
                           % (self.lnum,self.fname,self.lnum/dur))
            return affyrec

        def __iter__(self):
            return self

    def __init__(self, dbname='dbsnp', tablename='affy500k', annoflist = [],
                 addnullpk = 1):
        """
        """
        self.addnullpk = addnullpk # controls autoincrement primary key
        self.dbname = dbname
        self.tablename = tablename
        self.annoflist = annoflist
        try:
            # open annotation file as an iterable that returns a
            # decsved field list
            f = self.affyannoF(annoflist[0], addnullpk=addnullpk) 
        except BadHeaderError:
            LOG.error('The header for %s does not match our understanding ',
                   + 'of an Affymetrix 500k annotation file' % annoflist[0])
            sys.exit(1)
        fieldnames = f.colnames # freshly read and munged
        nfields = len(fieldnames)
        # sensible default for the lazy :-0
        fieldtypes = ['varchar(20)']*(nfields) 
        if self.addnullpk:
            # remember to insert a null here always
            fieldnames.insert(0,'ID') 
            # autoincrement pk if wanted
            fieldtypes.insert(0,'INTEGER PRIMARY KEY')
            nfields += 1
        self.fieldnames = fieldnames
        self.fieldtypes = fieldtypes
        self.createvars = ','.join(['%s %s' % (self.fieldnames[x],self.fieldtypes[x]) for x in range(nfields)])
        varplaces = ','.join(['?']*nfields) # sqlite uses qmarks
        self.insertsql = ('INSERT INTO %s (%s) values (%s)' 
                          % (self.tablename,','.join(self.fieldnames),varplaces))
        if debug:
            LOG.debug('AffyAnnoF initialized')
            LOG.debug('Got fieldnames = %s' % fieldnames)
            LOG.debug('createvars = %s' % self.createvars)
            LOG.debug('insertsql = %s' % self.insertsql)

           
    def buildSQLite(self):
        """create and fill a sqlite db from the mapping annotation files
        downloaded from
        http://www.affymetrix.com/support/technical/byproduct.affx?product=500k
        """
        from pysqlite2 import dbapi2 as sqlite
        # Create a connection to the database file "mydb":
        con = sqlite.connect(self.dbname)
        # Get a Cursor object that operates in the context of Connection con:
        cur = con.cursor()
        try:
            sql = 'drop table %s' % self.tablename
            cur.execute(sql)
        except:
            if debug:
                LOG.debug('no %s to drop this time!' % self.tablename)
            else:
                pass
        sql = 'create table %s (%s);' % (self.tablename,self.createvars)
        cur.execute(sql)
        for fname in self.annoflist:
            f = self.affyannoF(fname, addnullpk=self.addnullpk)
            LOG.debug('addnullpk: %s' % str(f.addnullpk))
            # just pass a file iterator!@
            cur.executemany(self.insertsql, f) 
            con.commit() # only at end of each file to speed things up
        cur.close()
        con.close()

    def buildMySQL(self,host, userid, password):
        """create and fill a mysql db from the mapping annotation files
        downloaded from
        http://www.affymetrix.com/support/technical/byproduct.affx?product=500k
        """
        import MySQLdb
        ntoinsert = 500
        # Create a connection to the database file "mydb"
        genome = MySQLdb.Connect('godzilla', 'refseq', 'Genbank')
        curs = genome.cursor() # use default cursor
        curs.execute('use %s' % self.dbname)
        varss = '''id INT UNSIGNED NOT NULL AUTO_INCREMENT, affyid varchar(15),
        probesetid varchar(15), rs varchar(15), chrom char(2), offset int(12), 
        strand char(1), index rsindex (rs), index chromrs (chrom, rs), 
        primary key (id) '''
        ivar = "id,affyid,probesetid,rs,chrom,offset,strand"
        qvall = ['%s']*7
        qvals = ','.join(qvall)
        insertsql = ("""insert into %s (%s) values (%s)""" 
                     % (self.tablename,ivar,qvals))
        # start again each time
        sql = ('DROP TABLE IF EXISTS %s.%s' 
               % (self.dbname,self.tablename)) 
        curs.execute(sql)
        sql = ('CREATE TABLE %s.%s (%s)' 
               % (self.dbname,self.tablename,varss))
        curs.execute(sql)
        for fname in self.annoflist:
            f = affyannoF(fname)
            done = 0
            insertme = []
            while not done:
                while (len(insertme) < ntoinsert) and not done:
                    try:
                        ll = f.next()
                        insertme.append(ll)
                    except StopIteration:
                        done = 1
                if len(insertme) > 0:
                    # just pass the list
                    curs.executemany(insertsql,insertme) 
                    insertme = []
            # only at end of each file to speed things up
            genome.commit() 
        curs.close()
        genome.close()

    def build(self, useMySQL = 0, host = '', userid = '', password = ''):
        if useMySQL:
            self.buildMySQL(host=host, userid=userid, password=password)
        else:
            self.buildSQLite()
            

def test(annofiles=[], pk=0, useMySQL=0):
    """exercise some basic functions during development"""
    affy = affyanno(dbname = 'AFFY500.sqlite.db', tablename='snpinfo', 
                    annoflist=annofiles, addnullpk = pk)
    affy.build(useMySQL = useMySQL)


if __name__ == "__main__":
    if len(sys.argv) > 1 :
        fspec = sys.argv[1]
        annofiles = glob.glob(fspec)
        if len(annofiles) > 0:
            annofiles.sort()
            test(annofiles, pk=0)
    else:
        print ('Please supply a wildcard in quotes for your affy mapping '
               + 'csv files - %s matches nothing' % fspec)
