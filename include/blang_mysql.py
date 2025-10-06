"""MySQL utility functions & package imports"""

# MySQL
# import mysql.connector as my
from blang import Die, rx, tq
import sqlalchemy as sa
# from pymysql.constants import CLIENT  # To enable multi-statement queries
import re
import sys

# Initialize
# By default: don't show MySQL queries as they are run (not "loud")
blang_superloudmysql = 0

# MySQL functions
def Connect(database = "blang", server = ""):
    global blang_mysql_connection
    
    if (server == ""):
        from sqlalchemy.engine.url import URL

        myDB = URL.create(drivername='mysql+pymysql', host='[mysql server]', # Insert MySQL server name
            database=database,
            query={
                'autocommit': '1',      # For InnoDB
                'local_infile': '1',    # For LOAD DATA LOCAL INFILE
                'read_default_file': '/home/blang1/.my.cnf'
            }
        )
        # blang_mysql_connection = sa.create_engine(url=myDB)
        blang_mysql_engine = sa.create_engine(url=myDB)
        # Multi-statement queries:
        # These are very difficult to implement with sqlalchemy. I haven't figured out how to retrieve results from anything but the first statement.
        # Meanwhile, getting warnings only seems possible by running SHOW WARNINGS, which presumably only gets the warnings for the last statement.
        # Disabling multi-statement queries (I very rarely use them).
        # blang_mysql_engine = sa.create_engine(url=myDB, connect_args={'client_flag': CLIENT.MULTI_STATEMENTS})    # To enable multi-statement queries
        # blang_mysql_engine = sa.create_engine(url=myDB, connect_args={'client_flag': CLIENT.MULTI_STATEMENTS}, echo=True)    # To enable multi-statement queries and print verbose output
        # blang_mysql_engine = sa.create_engine(url=myDB, connect_args={'local_infile': 1})

        blang_mysql_connection = blang_mysql_engine.connect()

        # # Always show warnings
        # Doesn't work with pymysql
        # def after_execute(conn, cursor, statement, parameters, context):
        #     warnings = blang_mysql_engine.dialect.serializer.get_warnings(result)
        #     for warning in warnings:
        #         print(warning)  # Process or log the warning
        # sa.event.listen(blang_mysql_engine, 'after_execute', after_execute)

    else:
        raise Exception(f"\n\nError: Unhandled server '{server}'\n")

# # Define custom derived class to be able to do len(query) on CursorResult objects
# This is going to be too inefficient since the CursorResult object gets copied. Not doing this just for adding len(query).
# class CursorResultWithLen(sa.engine.cursor.CursorResult):
# 
#     # def __init__(self, *args, **kwargs):
#     #     super().__init__(*args, **kwargs)
# 
#     # With composition
#     def __init__(self, query):
#         self.query = query
#         self.rowcount = self.query.rowcount
#         vars(self).update(vars(query))
# 
#     def __len__(self):
#         return self.query.rowcount
#     
#     def fetchall(self):
#         return self.query.fetchall
#     
#     def fetchone(self):
#         return self.query.fetchone
    
def Query(query, loud = None):
    global blang_mysql_connection
    global blang_superloudmysql
    
    # "loud" is a sticky parameter: once set for one query, all queries will be printed as well as run.
    # Explicitly use loud=0 on a query to disable this again.
    if loud == 1:
        blang_superloudmysql = 1
    if loud == 0:
        blang_superloudmysql = 0
    
    # Print query if in "loud" mode
    if blang_superloudmysql == 1:
        print("\n" + query)

    # c = blang_mysql_connection.cursor(buffered=True)
    # c = blang_mysql_connection.execute(sa.text(query))
    try:
        c = blang_mysql_connection.execute(sa.text(query))
    except:
        raise Exception(f"\n\nError: Query failed for query:\n\n{query}\n")
    
    warnings = blang_mysql_connection.execute(sa.text("SHOW WARNINGS"))
    if warnings.rowcount > 0:
        for warning in warnings:
            print(warning, file=sys.stderr)
        # Die on warnings
        # raise Exception(f"\n\nError: Query produced warnings for query:\n\n{query}\n")
        # Warn only
        print(f"\nWarning: Query produced warnings for query:\n\n{query}\n\n", file=sys.stderr)
    # warnings.close()

    # # Use custom CursorResult class with __len__
    # c = CursorResultWithLen(c)

    # Add query string to cursor object
    c.blang_query_string = query

    return c
    
def Numrows(query):
    return(query.rowcount)

def FetchList(query):
    """Fetch MySQL rows as a list (either of single values, or of tuples)"""
    a = query.fetchall()
    # If this is a single column:
    if (Numrows(query) > 0) and (max([len(x) for x in a]) == 1):
        # List of values
        return [x[0] for x in a]
    else:
        # List of tuples
        return a

def FetchSet(query):
    """Fetch MySQL rows as a set (either of single values, or of tuples)"""
    a = query.fetchall()
    # If this is a single column:
    if (Numrows(query) > 0) and (max([len(x) for x in a]) == 1):
        # Set of values
        return set([x[0] for x in a])
    else:
        # Set of tuples
        return set(a)

def FetchMap(query):
    """Fetch MySQL rows as a key-value mapping dictionary (must select at least two columns)"""
    a = query.fetchall()
    # Check if this is two columns
    if (Numrows(query) > 0) and min([len(x) for x in a]) == 2 and max([len(x) for x in a]) == 2:
        # Construct dictionary
        return dict([x for x in a])
    else:
        raise Exception(f"\n\nError: Expected 2 columns, but got between '{min([len(x) for x in a])}' and '{max([len(x) for x in a])}' for query:\n\n{query.blang_query_string}\n")

def FetchPanda(query):
    """Fetch MySQL rows as a pandas data frame using pd.read_sql()"""
    import pandas as pd
    global blang_mysql_connection
    try:
        df = pd.read_sql(query, blang_mysql_connection)
    except:
        raise Exception(f"\n\nError: Panda query failed for query:\n\n{query.blang_query_string}\n")
    return df
# Add alias
Panda = FetchPanda

def InsertPanda(panda, table, schema = None):
    """Insert Pandas data frame into SQL table"""
    global blang_mysql_connection
    
    # Parse schema name if defined in table name
    if schema is None:
        m = rx(r"^([^\.]+)\.([^\.]+)$", table)
        if m:
            if len(m) != 2:
                Die(f"Error: Couldn't parse schema from table name '{table}'")
            else:
                schema = m[0]
                table = m[1]
    
    # https://pandas.pydata.org/docs/reference/api/pandas.DataFrame.to_sql.html
    # DataFrame.to_sql(name, con, schema=None, if_exists='fail', index=True, index_label=None, chunksize=None, dtype=None, method=None)[source]
    return panda.to_sql(table, blang_mysql_connection, schema=schema, if_exists='append', index=False)
# Add alias
Insert = InsertPanda

def Fetch(query):
    return tq(query, total=Numrows(query))
# Alias
tqq = Fetch

def FetchRow(query):
    # Fetch a single MySQL row
    # Not really needed, can simply do e.g.:
    # query = Query(f"SELECT seq FROM …")
    # for (seq,) in tq(query, total=Numrows(query)):
    
    res = query.fetchone()

    # Convert tuple to str if there's only one field being returned (to avoid having to do e.g. "(name,) = Fetch(query)")
    if (len(res) == 1):
        (res,) = res

    # Return
    return res

def FetchDict(query):
    """Fetch a single MySQL row as a key-value dictionary (using the column names as keys)"""
    row = query.fetchone()
    # names = query.description
    names = query.keys()  # sqlalchemy
    # Get first element (the name) from each column's description
    names = [name[0] for name in names]
    # Construct dictionary (from a list of two-ples)
    res = dict([(names[i], row[i]) for i in range(len(names))])
    return res

def FetchOne(query):

    # Buffered (need this anyway to do Numrows(query))
    if query.rowcount == 1:
        res = query.fetchone()

        # Convert tuple to str if there's only one field being returned (to avoid having to do e.g. "(name,) = FetchOne(query)")
        if (len(res) == 1):
            (res,) = res
            # Convert str to int if numeric
            if type(res) == str:
                if rx(r'^\d+$', res):
                    res = int(res)

        # Return
        return res
    else:
        # sys.exit(f"\nError: Numrows {query.rowcount} instead of 1 for query:\n{query.statement}\n\nWarnings:\n{query.fetchwarnings()}\n");
        # sys.exit(f"\n\nError: Numrows {query.rowcount} instead of 1 (query: \"{query.statement}\")");
        # raise Exception(f"\nError: Numrows {query.rowcount} instead of 1 (query: \"{query.statement}\")");
        # print("\n\n")
        # raise Exception(f"\n\n\nError: Numrows {query.rowcount} instead of 1 for query:\n\n{query.statement}\n")
        # breakpoint()
        raise Exception(f"\n\nError: Numrows {query.rowcount} instead of 1 for query:\n\n{query.blang_query_string}\n")

    # # Unbuffered (where query.rowcount isn't available until all results are fetched):
    # # Fetch a row
    # res = query.fetchone()
    # # Check if it's empty
    # if res is None:
    #     # sys.exit(f"\nError: Numrows 0 instead of 1 for query:\n{query.statement}\n\nWarnings:\n{query.fetchwarnings()}\n");
    #     sys.exit(f"\nError: Numrows 0 instead of 1 (query: {query.statement})\n");
    # else:
    #     res2 = query.fetchone()
    #     if res2 is None:
    #         if (len(res) == 1):
    #             # Convert tuple to str if there's only one field being returned (to avoid having to do e.g. "(name,) = FetchOne(query)")
    #             (res,) = res
    #             return(res)
    #         else:
    #             return(res)
    #     else:
    #         # sys.exit(f"\nError: Numrows >1 instead of 1 for query:\n{query.statement}\n\nWarnings:\n{query.fetchwarnings()}\n");
    #         sys.exit(f"\nError: Numrows >1 instead of 1 (query: {query.statement})\n");
        
def Nextset(query):
    """Skip to the next result set returned by a multi-statement query"""
    query.nextset()

def Esc(s):
    """Escape a string for SQL"""

    # Escape \ to \\
    s = s.replace("\\", "\\\\")

    # Escape ' to \'
    s = s.replace("'", "\\'")

    # Escape : to \: (sqlalchemy-specific)
    s = s.replace(":", "\\:")

    return(s)

# import re
# def Esc(s):
#     """Escape string using mysql.connector function"""
#     global blang_mysql_connection
#     return blang_mysql_connection.converter.escape(s)

def Load(file, table, fields, silent=False):
    """Load TSV file into SQL table"""

    if silent == False: print(f"\nLoading '{file}' into table '{table}'...")
    Query(f"ALTER TABLE {table} DISABLE KEYS")
    Query(f"LOAD DATA LOCAL INFILE '{file}' INTO TABLE {table} ({', '.join(fields)}) SET id=NULL")

    if silent == False: print(f"Enabling keys...")
    Query(f"ALTER TABLE {table} ENABLE KEYS")

    if silent == False: print("Done!")

def Clear(table):
    global blang_mysql_connection
    print(f"\nClearing table '{table}'...\n");
    # c = blang_mysql_connection.cursor(buffered=True)
    # c.execute(f"TRUNCATE {table}")
    blang_mysql_connection.execute(sa.text(f"TRUNCATE {table}"))
    # print(f"\nCleared table '{table}'");
    # print("Done!");
        
def Optimize(table):
    global blang_mysql_connection
    # c = blang_mysql_connection.cursor(buffered=True)
    print(f"\nOptimizing table '{table}'...");
    # c.execute(f"ANALYZE TABLE {table}")
    # c.execute(f"OPTIMIZE TABLE {table}")
    blang_mysql_connection.execute(sa.text(f"ANALYZE TABLE {table}"))
    blang_mysql_connection.execute(sa.text(f"OPTIMIZE TABLE {table}"))
    print("Done!\n");
        
def FetchAll(query):
    return query.fetchall()
    
# Connect to MySQL
Connect()
