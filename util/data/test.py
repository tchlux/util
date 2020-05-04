
from util.data.data import Data, flatten

# TODO:  Test for file that has a type digression happen on a line
#        with empty strings in it. This caused a crash [2019-12-17].
# TODO:  Test data with no names getting other data added in place,
#        should overwrite the names!


# Some tests for the data
def test_data():
    import os
    import numpy as np

    print("Testing Data...", end=" ")

    # ----------------------------------------------------------------
    a = Data()

    # Verify append
    a.append([1,"a"])
    a.append([2,"b"])
    a.append([3,"c"])

    # Verify add column (with edge case, none type)
    a.add_column([None,None,None])
    assert(a.types[-1] == type(None))

    # Reassign the missing value column to floats, verify type update
    a[a.names[-1]] = [1.2, 3.0, 2.4]
    assert(a.types[-1] == float)
    # WARNING: Need to check doing (<int>, <str>) assignment

    # Verify automatic type generation AND automatic type-casting
    a.append([1,2,3])
    assert(tuple(a[-1]) == tuple([1,"2",3.0]))
    assert(tuple(map(type,a[-1])) == (int, str, float))

    # Add some missing values.
    a.append([-1,None,None])

    #  Done modifying "a", rest of tests leave it constant.
    # ----------------------------------------------------------------

    # Verify that "index" works correctly.
    assert( a.index([1,"2",3.0]) == 3 )
    assert( a.index([-1,None,None]) == 4 )

    # Verify that slicing works on descirptors.
    assert(a.names[:-1] == ['0','1'])

    # Verify that slicing works on descriptors of views.
    b = a[:,1:]
    assert(b.names[1:] == ['2'])

    # Verify in-place addition of rows (different length data with same column names)
    b = a[:-2].copy()
    c = a[1:-2]
    b += c
    assert(tuple(b["0"]) == tuple([1,2,3,2,3]))

    # Verify in-place addition of new columns (same length data with new columns)
    b = a[:,:-1].copy()
    c = a[:,-1:].copy()
    c["3"] = range(len(c))
    b += c
    assert(tuple(b[-2]) == (1,"2",3.0,3))
    assert(tuple(b[-1]) == (-1,None,None,4))

    # Verify in-place addition of new columns (same length data with fewer columns)
    b = a[:,:-1].copy()
    c = a[:,:].copy()
    c["3"] = range(len(c))
    c += b
    assert(tuple(c[-2]) == (1,"2",None,None))
    assert(tuple(c[:,-1]) == tuple(range(5)) + (None,)*5)

    # Verify the addition of in place addition of new columns AND rows.
    b = a[:]
    b['3'] = range(len(b))
    c = a[:]
    c += b
    assert(tuple(c['3']) == tuple([None]*len(b) + list(range(len(c)-len(b)))))

    # Verify in-place addition of a data set with *some* duplicate
    # column names, but the same number of rows (adds columns).
    b = a[:,1:].copy()
    c = a[:,:-1].copy()
    c.names = range(3,3+c.shape[1])
    b += c
    assert(b.shape == (5,4))
    assert(tuple(b.names) == tuple(map(str,range(1,5))))
    assert(len(list(b[:,-1] == c[:,-1])) == b.shape[0])

    # Verify slicing-based singleton value assignment
    b = a[:-2].copy()
    b[:,0] = 1
    assert(tuple(b["0"]) == tuple([1,1,1]))

    # Verify slicing-based multiple value assignment
    b = a[:-2].copy()
    b[:,0] = [3,1,4]
    assert(tuple(b["0"]) == tuple([3,1,4]))

    # Verify double indexing with integers
    assert(a[0,1] == "a")

    # Verify double indexing with slices
    assert(tuple(a[::-1,1])[2:] == tuple(["c","b","a"]))

    # Verify standard index access
    assert(tuple(a[0]) == tuple([1,"a",1.2]))

    # Verify column-access and automatic column naming
    assert(tuple(a["0"])[:-2] == tuple([1,2,3]))

    # Verify slicing by index
    assert(tuple(a[:2][0]) == tuple(a[0]))

    # Verify slicing by names
    assert(tuple(a["0","2"][0]) == tuple([1,1.2]))

    # Verify that copies with "copy" are deep
    b = a.copy()
    b.retype([str,str,str])
    assert(a.types[0] != str)

    # Assert that the two forms of full-data slicing make copies.
    assert(id(a[:].data) != id(a.data))
    assert(id(a[:,:].data) != id(a.data))

    # Verify that copies with slicing are deep and that retype works.
    b = a[:]
    b.retype([str,str,str])
    assert(a.types[0] != str)

    # Verify that in-place addition of data with same names in
    # different order is handled correctly.
    b = a[:,:]
    c = a[:,['1','0','2']]
    b += c
    assert(tuple(b[:,1]) == 2*tuple(a[:,1]))

    # Verify that multiple transformations on data yield correct view.
    b = a[:]
    b += a
    b += a
    c = b[::2]
    c = c[['1', '2']].unique()
    assert(c.shape[0] == len(a))
    assert(tuple(c.names) == tuple(b.names[1:]))

    # Verify that sorting a data object and a view object both work.
    b = a[:]
    b.sort()
    assert(tuple(b[:,0]) == tuple(sorted(a[:,0])))
    b = a[:-1, [0,1]]
    b.sort()
    assert(tuple(b[:,1]) == tuple(sorted(a[:-1,1])))
    # Assert that "a" (the object being viewed) is unmodified.
    assert(tuple(a['1']) == ('a','b','c','2',None))

    # Verify the stack process.
    b = a[:]
    c = a[:]
    c['2'] += 1.0
    b += c
    b.stack('2')
    assert(tuple(b[0,-1]) == (a[0,-1], a[0,-1]+1.0))

    # Testing the "unstack" operation.
    b = a[:]
    b += a
    b.pop(-1)
    b.stack('2')
    # Convert the "stacked" column into generator ogbjects.
    b['2'] = ((v for v in row) for row in b['2'])
    b[-1,'2'] = None
    # Unstack, and verify that the expected output happens.
    b.unstack('2')
    assert(tuple(b['1']) == ('a','a','b','b','c','c','2','2',None))

    # Verify that reorder works.
    b = a[:]
    b.reorder(["2","0"])
    assert(tuple(b.names) == ("2","0","1"))
    assert(tuple(b.types) == (float,int,str))
    for i in range(len(b)):
        for n in b.names:
            assert(a[i,n] == b[i,n])

    # Verify the ability to construct real-valued arrays (and go back)
    numeric = a.to_matrix()
    c = Data(map(numeric.from_real, numeric))
    assert(tuple(a[:-1]["1"]) == tuple(c["1"]))

    # Verify conversion into numpy structured array (without "None")
    b = a.to_struct()
    assert(type(b) == np.ndarray)
    assert(type(b.dtype.names) == tuple)
    assert(tuple(b['1']) == tuple(a['1'][:-1]))

    # Verify column item retrieval and assignment
    assert(a["0"][0] == a[0][0])
    a["0"][0] = -10
    a["0"][0] = 1

    # Verify total column assignment
    first_col = list(a["0"])
    new_col = list(map(str,a["0"]))
    a["0"] = new_col
    assert(tuple(a["0"]) == tuple(new_col))
    a["0"] = first_col

    # Verify new column assignment
    b = a[:]
    first_col = list(b["0"])
    new_col = list(map(str,b["0"]))
    b["0-str"] = new_col
    assert(tuple(map(str,b["0"])) == tuple(b["0-str"]))

    # Verify copying a data object that has a 'None' in the first row
    b = Data(names=["0","1"], types=[int,str])
    b.append([0,None])
    b.append([1,'b'])
    c = b[:]
    assert(tuple(b[0]) == tuple(c[0]))

    # Verify accessing multiple rows by list-access
    b = a[[0,2,3]]
    assert(tuple(b["0"]) == tuple([1,3,1]))

    # Verify load and save of a csv file
    a.save("a-test.csv")
    b = Data().load("a-test.csv")
    assert(tuple(a["0"]) == tuple(b["0"]))
    os.remove("a-test.csv")

    # Verify the writing of a CSV with quoted content and correct read.
    b = a[:]
    b[-1,1] = "string with comma, we'll see"
    b.save("a-test.csv")
    c = Data.load("a-test.csv")
    assert(tuple(c[-1]) == tuple(b[-1]))
    os.remove("a-test.csv")

    # Verify load and save of a pkl file
    a.save("a-test.pkl")
    b = Data.load("a-test.pkl")
    assert(tuple(a["0"]) == tuple(b["0"]))
    os.remove("a-test.pkl")

    # Verify load and save of a gzipped dill file
    a.save("a-test.dill.gz")
    b = Data().load("a-test.dill.gz")
    assert(tuple(a["0"]) == tuple(b["0"]))
    os.remove("a-test.dill.gz")

    # Verify load and save of a gzipped csv file
    a.save("a-test.csv.gz")
    b = Data().load("a-test.csv.gz")
    assert(tuple(a["0"]) == tuple(b["0"]))
    os.remove("a-test.csv.gz")

    # Verify the behavior of the 'unpack' and 'pack' functions.
    b = a[:]
    b['iter'] = [None] + [list(range(i)) for i in range(len(b)-1)]
    _ = b.shape[1]
    b.unpack("iter")
    assert((b.shape[1] - _) == (len(b) - 3))
    b.pack("iter")
    assert(b.shape[1] == _)
    assert(b[0,-1] is None)
    assert(b[1,-1] is None)

    # Check the 'flatten' function to make sure it correctly flattens
    # structures with nested, differential depth, and different types.
    complex_sequence = [(l for l in ((1,2,3)+tuple("abc"))),
                        "def", "ghi", ("jk",tuple("lmn"))]
    simplified_sequence = (1,2,3)+tuple('abcdefghijklmn')
    assert(tuple(flatten(complex_sequence)) == simplified_sequence)

    # TODO: Verify load of a large file (sample only)
    # TODO: Verify load of a large compressed file (sample only)

    # Verify column equivalence / equality checker by element
    b = a[a["0"] == 1]
    assert(tuple(b["0"]) == tuple([1,1]))

    # Verify one of the inequality operators
    b = a[a["0"] < 2]
    assert(tuple(b["0"]) == tuple([1,1,-1]))

    # Verify column membership checker by set
    b = a[a["1"] == {'a','b','2'}]
    assert(tuple(b["1"]) == tuple(['a','b','2']))

    # Verify column index can be an integer or a string.
    for i in range(len(a)):
        for j,name in enumerate(a.names):
            assert(a[i,j] == a[i,name])

    # Verify that the "in" operator works on rows.
    assert(1 in a[0])
    assert('a' in a[0])
    assert(1.2 in a[0])

    # Verify that "in" works on Data and View.
    assert(not ([1, "a"] in a))
    assert(not ([3, "a"] in a[:,:2]))
    assert([1, "a"] in a[:,:2])

    # Verify that "equality" works for rows.
    assert(a[0] == (1,'a',1.2))
    assert(not all(a[0] != a[0][:]))

    # Verify that "equality" works for Data objects.
    assert(tuple(a == [1,'a',1.2]) == (0,))

    # Verify that "equality" fails when the equality operator for
    # elements of a row is NotImplemented.
    b = a[:]
    assert( len(list(b == (i for i in (1,2,3,4)))) == 0 )
    assert( tuple(b[0] == (i for i in (1,2,3,4))) ==
            (NotImplemented, NotImplemented, NotImplemented) )


    # Verify that "addition" with sequences works correctly (left and right).
    b = a[:]
    b[:,1] = map(lambda v: (ord(v) if (v != None) else None), b[:,1])
    assert(sum(b[0] + [-1,-97,-1.2]) == 0)
    assert(sum([-1,-97,-1.2] + b[0]) == 0)
    assert(tuple(a[0] + a[0]) == (2,'aa',2.4))
    # Verify that "subtraction" with sequences works correctly (left and right).
    assert(sum(b[0] - [1,97,1.2]) == 0)
    assert(sum([1,97,1.2] - b[0]) == 0)

    # WARNING: Not individually verifying *all* comparison operators.

    # Verify generator-based index assignment.
    b = a[:]
    b[(i for i in range(0,len(b),2)), "1"] = "new"
    assert(b[0,"1"] == b[2,"1"] == b[4,"1"] == "new")

    # Verify list-based generator assignment.
    b = a[:]
    b[[i for i in range(0,len(b),2)], "1"] = str("test")
    assert(b[0,"1"] == b[2,"1"] == b[4,"1"] == "test")

    # Verify that assignment from a generator works.
    b = a[:]
    b['3'] = map(str, b['2'])
    assert(b[0,'3'] == str(b[0,'2']))
    assert(b[-1,'3'] == str(b[-1,'2']))

    # Verify assignment of a non-iterable being broadcast automatically.
    b['4'] = None
    assert(tuple(b[-2]) == (1,'2',3.0,'3.0',None))

    # Vefify addition of a column to empty data works.
    b = Data()
    b.add_column([1,2,3])
    assert(tuple(b[:,0]) == (1,2,3))

    # Verify the generation of a random k-fold cross validation.
    b = Data()
    for i in range(37):
        b.append([i])
    # Collect the list of all the "testing" rows seen
    all_test = []
    for (train, test) in b.k_fold(k=11, seed=1):
        all_test += list(test["0"])
    assert(sorted(all_test) == list(range(len(b))))

    # Verify the identification of unique rows and collection process.
    b = a[:]
    b += a
    cols = ['0','1']
    b = b[cols].unique().copy().collect(b)
    assert(tuple(b[0,-1]) == tuple(2*[a[0,-1]]))

    # Verify slicing data down to one column.
    assert(tuple(a[:,-1]) == (1.2,3.0,2.4,3.0,None))

    # Test slicing to access columns and rows.
    b = a[:,:-1]
    c = b.unique()
    assert(len(b) == len(c))

    # Test the access of a single column by name in a view.
    b = a[:,:-1]
    assert( tuple(b[:,"0"]) == (1,2,3,1,-1) )

    # Verify that the printed data set handles edge cases (one more
    # row or equal rows to number desired) correctly.
    b = a[:]
    b = b + b
    b.append(b[0])
    sb = str(b)
    # Test the correctness of truncated rows and columns.
    b.max_display = 2
    b_printout = '''
========================
Size: (11 x 3)          
                        
 0   | ... | 2          
 int | ... | float      
------------------      
 1   | ... | 1.2        
 ...   ...   ...        
 1   | ... | 1.2        
                        
 missing 4 of 33 entries
      at rows: [5, 10]  
   at columns: [1, 2]   
========================
'''
    assert( str(b) == b_printout)
    # Test the correctness of truncated rows only.
    b.max_display = 4
    b_printout = '''
========================
Size: (11 x 3)          
                        
 0   | 1    | 2         
 int | str  | float     
-------------------     
 1   | "a"  | 1.2       
 2   | "b"  | 3.0       
 ...   ...    ...       
 -1  | None | None      
 1   | "a"  | 1.2       
                        
 missing 4 of 33 entries
      at rows: [5, 10]  
   at columns: [1, 2]   
========================
'''
    assert( str(b) == b_printout)
    # Test the correctness of exact match (no truncation).
    b.max_display = 11
    b_printout = '''
========================
Size: (11 x 3)          
                        
 0   | 1    | 2         
 int | str  | float     
-------------------     
 1   | "a"  | 1.2       
 2   | "b"  | 3.0       
 3   | "c"  | 2.4       
 1   | "2"  | 3.0       
 -1  | None | None      
 1   | "a"  | 1.2       
 2   | "b"  | 3.0       
 3   | "c"  | 2.4       
 1   | "2"  | 3.0       
 -1  | None | None      
 1   | "a"  | 1.2       
                        
 missing 4 of 33 entries
      at rows: [5, 10]  
   at columns: [1, 2]   
========================
'''
    assert( str(b) == b_printout)

    # Test the 'fill' method.
    b = a[:]
    b[-2,1] = 'd'
    b.append( b[0] )
    assert(str(b.fill(b[-2]).data) == str([-1,'a',1.7999999999999998]))
    assert(str(a.fill(a[-1][:]).data) == str([-1,'2',2.1]))

    # Test cases for using "fill" with weights.
    b.pop(-2)
    b["3"] = b["2"]
    b["2"] = b["1"]
    b[-1,-2] = 'b'
    b[-1,-1] = None
    condense = lambda v: str(v)[:4]
    assert(tuple(map(condense, b.fill(b[-1], weights=[1., 2., 3., 4.]).data))
           == ('1', 'a', 'b', '2.50'))
    b[-1,-1] = None
    assert(tuple(map(condense, b.fill(b[-1], weights=[100., 10., 1., 1000.]).data))
           == ('1', 'a', 'b', '1.21'))

    # Verify that the names and types attributes are the correct type.
    assert("Descriptor" in str(type(a.types)))
    assert("Descriptor" in str(type(a.names)))

    # ----------------------------------------------------------------
    #     Testing expected exceptions.

    # Verify that indexing an element that doesn't exist fails.
    try: a.index([1,2,3])
    except ValueError: pass
    else: assert(False)

    # Verify that incorrect names are detected in "Data.effect".
    b = a[:]
    try: b.effect("bad name")
    except Data.BadIndex: pass
    else:  assert(False)

    # Verify that a data object refuses to be added to itself in place.
    b = a[:]
    try: b += b
    except Data.Unsupported: pass
    else:  assert(False)

    # Try assigning with a generator that is too short.
    try:   a['3'] = (i for i in range(len(a)-1))
    except Data.BadData: pass
    else:  assert(False)

    # Try assigning with a generator that is too long.
    try:   a['3'] = (i for i in range(len(a)+1))
    except Data.BadData: pass
    else:  assert(False)

    # Try providing a generator that does not give strictly integers
    try:   a[(v for v in ('a',1,'b'))]
    except Data.BadIndex: pass
    else: assert(False)

    # Try adding a too-short row
    b = a[:]
    try:   b.append([9,"z"])
    except Data.BadElement: pass
    else:  assert(False)

    # Try a too-short column
    b = a[:]
    try:   b.add_column([1])
    except Data.BadData: pass
    else:  assert(False)

    # Try a too-long column
    b = a[:]
    try:   b.add_column(map(float,range(1000)))
    except Data.BadData: pass
    else:  assert(False)

    # Try adding an empty column
    b = a[:]
    try:   b.add_column([])
    except Data.BadData: pass
    else:  assert(False)

    # Try adding a name that is not a string
    b = a[:]
    try:   b.add_column([1,2,3,4,5], name=1)
    except Data.BadSpecifiedName: pass
    else:  assert(False)

    # Try adding a column with multiple types
    b = a[:]
    try:   b.add_column([1,2,3,4,5.0])
    except Data.BadValue: pass
    else:  assert(False)

    # Try a mismatched-type slice assignment
    b = a[:]
    try:   b[:,0] = [1.0, 1, 1, 1, 1]
    except Data.BadAssignment: pass
    else:  assert(False)

    # Try a mismatched-type column assignment operation
    try:   a["0"] = list(a["0"])[:-1] + ["0"]
    except Data.BadAssignment: pass
    else:  assert(False)

    # Try a too-short column assignment operation
    try:   a["0"] = list(map(str,a["0"]))[:-1]
    except Data.BadValue: pass
    else: assert(False)

    # Try a too-long column assignment operation
    try:   a["0"] = list(map(str,a["0"])) + [1]
    except Data.BadAssignment: pass
    else:  assert(False)

    # Test for an error when an added data set has duplicate columns.
    b = a[:]
    c = a[:]
    c.names[2] = '1'
    try:    b.collect(c)
    except: Data.BadData
    else:   assert(False)

    # Try a bad combination of names and types
    try:   Data(names=["a","b","c"], types=[str, str])
    except Data.NameTypeMismatch: pass
    else:  assert(False)

    # Try a bad value in "names"
    try:   Data(names=[1,2,3], types=["a","b","c"])
    except Data.BadSpecifiedName: pass
    else:  assert(False)

    # Try a bad value in "types"
    try:   Data(names=["a","b","c"], types=[1,2,3])
    except Data.BadSpecifiedType: pass
    else:  assert(False)

    # Try bad value
    try:   a.append(["a","b","c"])
    except Data.BadValue: pass
    else:  assert(False)

    # Try bad element
    try:   a.append([1,2])
    except Data.BadElement: pass
    else:  assert(False)

    # Try non-existing column index
    try:   a["hello"]
    except Data.UnknownName: pass
    else:  assert(False)

    # Attempt to set item on an empty data
    b = Data()
    try:   b[0] = 1
    except Data.Empty: pass
    else:  assert(False)

    # Try get item on an uninitialized data
    try:   Data()["a"]
    except Data.Empty: pass
    else:  assert(False)

    # Try retyping an uninitialized data
    try:   Data().retype([])
    except Data.Empty: pass
    else:  assert(False)

    # Try reordering an uninitialized data
    try:   Data().reorder([])
    except Data.Empty: pass
    else:  assert(False)

    # Try copying an uninitialized data
    try:   Data().copy()
    except Data.Empty: pass
    else:  assert(False)

    # Try popping from an uninitialized data
    try:   Data().pop()
    except Data.Empty: pass
    else:  assert(False)

    # Try saving an uninitialized data
    try:   Data().save("")
    except Data.Empty: pass
    else:  assert(False)

    # Done testing
    print("passed.")



# =============================================================
#      END OF Python 'Data' with named and typed columns     
# =============================================================


# Run test cases
if __name__ == "__main__":
    # Get the path to be watching for code coverage.
    import os
    path_to_watch = os.path.abspath(os.curdir)
    # Activate a code coverage module.
    import coverage
    cov = coverage.Coverage(source=[path_to_watch])
    cov.start()

    # Run the tests.
    test_data()

    # Save data from the coverage tracking (for report).
    cov.stop()
    cov.save()
    # Create a temporary directory for holding test results.
    from tempfile import TemporaryDirectory
    temp_dir = TemporaryDirectory()
    results_dir = temp_dir.name
    cov.html_report(directory=results_dir)
    # Open the results file.
    import webbrowser
    webbrowser.open("file://"+os.path.join(results_dir, "index.html"))
    # Wait for load, then delete the temporary directory.
    #   (user might decide to kill this process early).
    import time
    time.sleep(60*10)
    temp_dir.cleanup()
    del temp_dir
