'''
SBtab Tools
===========

These functions facilitate the use of SBtab. 
They can be used to create SBtab objects, by merging strings or read files, respectively.
'''

import tablib
import copy
import os.path
import sqlite3

from .SBtab import SBtabError, SBtabTable
from .SBtab import tablibIO


def oneOrMany(spreadsheet_file):
    '''
    Checks for multiple tables in a file and cuts them into separate tablib object.

    Parameters
    ----------
    spreadsheet_file : tablib object
        Tablib object of the whole SBtab table.
    '''
    sbtabs = []

    # Copy file, one for iteration, one for cutting
    sbtab_document = copy.deepcopy(spreadsheet_file)
    
    # Create new tablib object
    sbtab = tablib.Dataset()

    # Cutting sbtab_document, write tablib objects in list
    if len(spreadsheet_file) != 0:  # If file not empty
        for row in spreadsheet_file:
            if len(sbtab) == 0:  # If first line, append line w/o checking
                sbtab.rpush(sbtab_document.lpop())
            else:
                for i, entry in enumerate(row):
                    # If header row (!!), write to new tablib object and store the last one
                    if entry.startswith('!!'):
                        sbtabs.append(sbtab)
                        sbtab = tablib.Dataset()
                        sbtab.rpush(sbtab_document.lpop())
                        break
                    # If not header row, append line to tablib object
                    if len(row) == i + 1:
                        sbtab.rpush(sbtab_document.lpop())
        sbtabs.append(sbtab)

    # Return list of tablib objects
    return sbtabs


def openSBtab(filepath):
    '''
    Opens SBtab from file path.

    Parameters
    ----------
    filepath : str
        Path of the spreadsheet file.
    '''
    if not os.path.isfile(filepath):
        return None

    dataset = tablibIO.importSet(filepath)
    sbtab = SBtabTable(dataset, filepath)
    
    return sbtab


def openMultipleSBtab(filepath):
    '''
    Opens one or more SBtabTables from a single file.

    Parameters
    ----------
    filepath : str
        Path of the spreadsheet file.

    Returns
    ----------
    A list of SBtabTable objects.
    '''
    tablib_obj = tablibIO.importSet(filepath)
    datasets = oneOrMany(tablib_obj)
    return [SBtabTable(ds, filepath) for ds in datasets]


def openMultipleSBtabFromFile(f):
    '''
    Opens one or more SBtabTables from a file-like object.
    Assumes tab-separated.

    Parameters
    ----------
    f : file-like object.

    Returns
    ----------
    A list of SBtabTable objects.
    '''
    tablib_obj = tablibIO.haveTSV(f.read(), '\t')
    datasets = oneOrMany(tablib_obj)
    return [SBtabTable(ds, 'dummy.tsv') for ds in datasets]


def createDataset(header_row, columns, value_rows, filename):
    '''
    Creates an SBtab object by merging strings or list of strings.
    Takes a header row, main column row, and the value rows as lists of
    strings and returns an SBtab object.

    Parameters
    ----------
    header_row : str
        String of the header row.
    columns: list
        List of strings, names of the columns.
    value_rows : list
        List of lists containing the different rows of the table.
    '''
    # Initialize variables
    sbtab_temp = []
    sbtab_dataset = tablib.Dataset()
    header = header_row.split(' ')

    # Delete spaces in header, main column and data rows
    header = [x.strip(' ') for x in header]
    columns = [x.strip(' ') for x in columns]
    for row in value_rows:
        try:
            for entry in row:
                entry = entry.strip(' ')
        except:
            continue

    # Add header, main column and data rows to temporary list object
    sbtab_temp.append(header)
    sbtab_temp.append(columns)
    for row in value_rows:
        sbtab_temp.append(row)

    # Delete all empty entries at the end of the rows
    for row in sbtab_temp:
        if len(row) > 1:
            while not row[-1]:
                del row[-1]

    # Make all rows the same length
    longest = max([len(x) for x in sbtab_temp])
    for row in sbtab_temp:
        if len(row) < longest:
            for i in range(longest - len(row)):
                row.append('')
            sbtab_dataset.append(row)
        else:
            sbtab_dataset.append(row)

    # Create SBtab object from tablib dataset
    sbtab = SBtabTable(sbtab_dataset, filename)
    return sbtab

class SBtabDict(dict):
    
    def __init__(self, sbtab_list):
        """
            Arguments:
                sbtab_list - a list of SBtabTable objects
        """
        self.fpath = ''
        self.sbtab_list = sbtab_list
        for m in sbtab_list:
            self[m.table_name] = m

    @staticmethod
    def FromSBtab(fpath):
        spreadsheet_file = tablibIO.loadTSV(fpath, False)
        m = oneOrMany(spreadsheet_file)
        sbtab_list = [SBtabTable(dset, fpath) for dset in m]
        sbtab_dict = SBtabDict(sbtab_list)
        sbtab_dict.fpath = fpath
        return sbtab_dict

    def GetColumnFromTable(self, table_name, column_name):
        """
            Returns:
                a list of the values in the column called 'column_name'
                in the table 'table_name'
        """
        column_index = self[table_name].columns_dict['!' + column_name]
        rows = self[table_name].getRows()
        return [r[column_index] for r in rows]

    def GetColumnsFromTable(self, table_name, column_names):
        """
            Arguments:
                table_name   - the name of the table in the SBtab file (without '!!')
                column_names - a list of column names from which to get the data (without '!')
                
            Returns:
                a list of lists containing the values corresponding to the
                columns in 'column_names' in the table 'table_name'
        """
        try:
            idxs = [self[table_name].columns_dict['!' + c] for c in column_names]
        except KeyError as e:
            all_columns = ', '.join(self[table_name].columns_dict.keys())
            raise KeyError('Cannot find the column %s in table "%s" in file %s. '
                           'Columns are: %s'
                           % (e, table_name, self.fpath, all_columns))
        return [list(map(r.__getitem__, idxs)) for r in self[table_name].getRows()]
        
    def GetDictFromTable(self, table_name, key_column_name, value_column_name,
                         value_mapping=None):
        column_names = [key_column_name, value_column_name]
        key_val_list = self.GetColumnsFromTable(table_name, column_names)
        if value_mapping is not None:
            keys, vals = list(zip(*key_val_list))
            vals = list(map(value_mapping, vals))
            key_val_list = zip(keys, vals)
        return dict(key_val_list)
        
    def GetTableAttribute(self, table_name, attribute_name):
        """
            Arguments:
                table_name     - the name of the table in the SBtab file (without '!!')
                attribute_name - a string with the attribute name
                
            Returns:
                A string containing the value of the attribute in that table,
                or None if the attribute does not exist
        """
        try:
            return self[table_name].getCustomTableInformation(attribute_name)
        except SBtabError:
            return None
            
    def SBtab2SQL(self, comm, append=False):
        comm.execute("CREATE TABLE IF NOT EXISTS __tables__ (TableName TEXT, TableType TEXT, "
                     "header TEXT)")
        comm.execute("CREATE TABLE IF NOT EXISTS __columns__ (TableName TEXT, idx INT, ColumnName TEXT)")

        for m in self.sbtab_list:
            # get the names of the columns in the right order (i.e. so that
            # the corresponding column indices will be 0..n)
            
            columns = sorted(m.columns, key=m.columns_dict.get)
            columns = map(lambda c: str(c[1:]), columns)
            columns = [c for c in columns if c != '']
            
            rows = list(comm.execute("SELECT * FROM __tables__ "
                                     "WHERE TableName = '%s'" %  m.table_name))
            if len(rows) > 0:
                # if the table already exists, make sure that the metadata is
                # the same as in the SBtab.
                tname, ttype, theader = rows[0]
                assert ttype == m.table_type
                assert theader == m._getHeaderRow()
                
                # TODO: also assert that the columns are exactly the same as before
                
                if not append:
                    comm.execute("DROP TABLE %s" % m.table_name)
            else:
                # if the table doesn't already exist, add an entries for it 
                # in the __tables__ and __columns__
                comm.execute("INSERT INTO __tables__ VALUES(?,?,?)", 
                             [m.table_name, m.table_type, m._getHeaderRow()])
    
                for i, col in enumerate(columns):
                    comm.execute("INSERT INTO __columns__ VALUES(?,?,?)", 
                                 [m.table_name, i, col])
            
            col_text = ','.join(['\'%s\' TEXT' % col for col in columns])
            comm.execute("CREATE TABLE IF NOT EXISTS %s (%s)" % (m.table_name, col_text))

            # copy the data from the SBtab table into the relevant table in the 
            # database.
            ins_command = "INSERT INTO %s VALUES(%s)" % \
                          (m.table_name, ','.join(["?"]*len(columns)))
            for i, row in enumerate(m.getRows()):
                if len(row) > len(columns):
                    row = row[0:len(columns)]
                comm.execute(ins_command, row)
        
        comm.commit()

    @staticmethod
    def FromSQLite(fpath):
        """
            Read all tables from a SQL database into an SBtab object.
            This function assumed that the database has one table
            called __tables__ with the relevant header fields for SBtab
        """
        comm = sqlite3.connect(fpath)
        assert list(comm.execute("SELECT name FROM sqlite_master WHERE name='__tables__'")) != []
        table_names, table_types, headers = \
            zip(*comm.execute("SELECT TableName, TableType, header from __tables__"))
        
        sbtabs = []
        for table_name, header in zip(table_names, headers):

            columns = []
            for c in comm.execute("SELECT ColumnName from __columns__ WHERE "
                                  "TableName == '%s' ORDER BY idx" % table_name):
                columns.append(c[0])
            
            sbtab = tablib.Dataset()
            sbtab.rpush([header] + [''] * (len(columns)-1))
            sbtab.rpush(map(lambda s: '!' + s, columns))
            for row in comm.execute("SELECT * FROM '%s'" % table_name):
                sbtab.append(row)
            sbtabs.append(sbtab)

        sbtab_list = [SBtabTable(dset, fpath) for dset in sbtabs]
        sbtab_dict = SBtabDict(sbtab_list)
        sbtab_dict.fpath = fpath
        comm.close()
        return sbtab_dict
