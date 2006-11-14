import sys, glob, time
import logging

debug = True

logging.basicConfig(level=logging.DEBUG)
LOG = logging.getLogger("readAffySeqFasta")


class BadHeaderError(Exception):
    """catch bad headers
    """
    pass

      
class AffySeqFasta:
    """class to handle affy 500k mapping sequence fasta files and write as
    sqlite table.
    """

    class AffySeqIter:
        """Iterable for fasta file that returns value suitable for 
        sqlite executemany.

        """

        def __init__(self, fname='', addnullpk=0):
            self.colnames = ['probe_set_id', 'sequence']
            self.f = file(fname,'r')
            self.lnum = 0
            self.fname = fname
            self.started = time.time()
            # flag for autoincrement pk
            self.addnullpk = addnullpk 

        def extract_id(self, id_line):
            if id_line.startswith('>affx:'):
                return id_line[6:]
            raise BadHeaderError("near line %d" % self.lnum)
        
        def next(self):
            try:
                ## The sequence fasta file looks like this:
                ##   >affx:SNP_A-1780271
                ##   ggataaaattagagaaSacatttgttgaaccat
                ##   >affx:SNP_A-1780274
                ##   ggataaagaagagaagRattcctggcagcagta
                id_line = self.f.next().strip()
                seq_line = self.f.next().strip()
                self.lnum += 2
            except:  ## FIXME: bare except is bad, find Exception class
                raise StopIteration
            rec = []
            rec.append(self.extract_id(id_line))
            rec.append(seq_line)
            if self.addnullpk:
                rec.insert(0, None) # for autoincrement pk
            if self.lnum % 10000 == 0:
                if debug:
                    dur = time.time() - self.started
                    print ('reading line %d of %s = %5.2f recs/sec' 
                           % (self.lnum, self.fname, self.lnum/dur))
            return rec

        def __iter__(self):
            return self

    def __init__(self, dbname='dbsnp', tablename='seq', flist = [],
                 addnullpk = 1):
        """
        """
        self.addnullpk = addnullpk # controls autoincrement primary key
        self.dbname = dbname
        self.tablename = tablename
        self.flist = flist
        try:
            # open annotation file as an iterable that returns a
            # decsved field list
            f = self.AffySeqIter(flist[0], addnullpk=addnullpk) 
        except BadHeaderError:
            LOG.error('The header for %s does not match our understanding ',
                   + 'of an Affymetrix 500k annotation file' % annoflist[0])
            sys.exit(1)
        fieldnames = f.colnames
        nfields = len(fieldnames)
        # sensible default for the lazy :-0
        fieldtypes = ['varchar(20)'] * (nfields) 
        if self.addnullpk:
            # remember to insert a null here always
            fieldnames.insert(0, 'ID') 
            # autoincrement pk if wanted
            fieldtypes.insert(0, 'INTEGER PRIMARY KEY')
            nfields += 1
        self.fieldnames = fieldnames
        self.fieldtypes = fieldtypes
        self.createvars = ','.join(['%s %s' % (self.fieldnames[x], self.fieldtypes[x]) for x in range(nfields)])
        varplaces = ','.join(['?'] * nfields) # sqlite uses qmarks
        self.insertsql = ('INSERT INTO %s (%s) values (%s)' 
                          % (self.tablename, ','.join(self.fieldnames), varplaces))
        if debug:
            LOG.debug('createvars = %s' % self.createvars)
            LOG.debug('insertsql = %s' % self.insertsql)

           
    def buildSQLite(self):
        """create and fill a sqlite db from the sequense fasta file
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
        sql = 'create table %s (%s);' % (self.tablename, self.createvars)
        cur.execute(sql)
        for fname in self.flist:
            f = self.AffySeqIter(fname, addnullpk=self.addnullpk)
            # just pass a file iterator
            cur.executemany(self.insertsql, f) 
            con.commit() # only at end of each file to speed things up
        cur.close()
        con.close()

    def build(self):
        self.buildSQLite()
            

def test(files=[], pk=0):
    """exercise some basic functions during development"""
    affy = AffySeqFasta(dbname = 'AFFY500.sqlite.db', tablename='seq',
                        flist=files, addnullpk = pk)
    affy.build()


if __name__ == "__main__":
    fspec = 'Mapping250K*fasta'
    if len(sys.argv) > 1 :
        fspec = sys.argv[1]
    annofiles = glob.glob(fspec)
    if len(annofiles) > 0:
        annofiles.sort()
        test(annofiles, pk=0)
    else:
        print ('Please supply a wildcard in quotes for your affy seq '
               + 'fasta files - %s matches nothing' % fspec)
