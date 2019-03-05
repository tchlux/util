from util.data.data import Data

# Some tests for the data
def test_data():
    import os
    import numpy as np

    print("Testing Data...", end=" ")

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

    # Verify in-place addition
    b = a[:]
    c = a[1:]
    b += c
    assert(tuple(b["0"]) == tuple([1,2,3,2,3]))

    # Verify slicing-based singleton value assignment
    b = a[:]
    b[:,0] = 1
    assert(tuple(b["0"]) == tuple([1,1,1]))

    # Verify slicing-based multiple value assignment
    b = a[:]
    b[:,0] = [3,1,4]
    assert(tuple(b["0"]) == tuple([3,1,4]))

    # Verify double indexing with integers
    assert(a[0,1] == "a")

    # Verify double indexing with slices
    assert(tuple(a[::-1,1]["1"]) == tuple(["c","b","a"]))

    # Verify standard index access
    assert(tuple(a[0]) == tuple([1,"a",1.2]))

    # Verify column-access and automatic column naming
    assert(tuple(a["0"]) == tuple([1,2,3]))

    # Verify slicing by index
    assert(tuple(a[:1][0]) == tuple(a[0]))

    # Verify slicing by names
    assert(tuple(a["0","2"][0]) == tuple([1,1.2]))

    # Verify that copies with "copy" are deep
    b = a.copy()
    b.retype([str,str,str])
    assert(a.types[0] != str)

    # Verify that copies with slicing are deep and that retype works.
    b = a[:]
    b.retype([str,str,str])
    assert(a.types[0] != str)

    # Verify that reorder works.
    b = a[:]
    b.reorder(["2","0"])
    assert(tuple(b.names) == ("2","0","1"))
    assert(tuple(b.types) == (float,int,str))
    for i in range(len(b)):
        for n in b.names:
            assert(a[i,n] == b[i,n])

    # Verify automatic type generation AND automatic type-casting
    a.append([1,2,3])
    assert(tuple(a[-1]) == tuple([1,"2",3.0]))
    assert(tuple(map(type,a[-1])) == (int, str, float))

    # Verify caching of missing values
    a.append([-1,None,None])
    assert(len(a.missing) == 1)

    # Verify the ability to construct real-valued arrays (and go back)
    numeric = a.to_matrix()
    c = Data(map(numeric.real_to, numeric.data))
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

    # Verify that missing values are handled correctly by add_column
    b = a[:]
    b.add_column([1,2,3,4,None])
    assert(id(b[-1]) in b.missing)
    assert(len(b.missing) == 1)

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
    b = b[cols].unique().collect(b)
    assert(tuple(b[0,-1]) == tuple(2*[a[0,-1]]))

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

    # Attempt to index an empty data
    b = Data()
    try:   b[0]
    except Data.Empty: pass
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

    # Try accessing by row-index and having duplicate entries.
    try:   a[[0,0]]
    except Data.BadIndex: pass
    else: assert(False)

    # Try adding a too-short row
    try:
        b = a[:]
        b.append([9,"z"])
    except Data.BadElement: pass
    else: assert(False)

    # Try a too-short column
    try:
        b = a[:]
        b.add_column([1])
    except Data.BadData: pass
    else: assert(False)

    # Try a too-long column
    try:
        b = a[:]
        b.add_column(map(float,range(1000)))
    except Data.BadData: pass
    else: assert(False)

    # Try adding a name that is not a string
    try:
        b = a[:]
        b.add_column([1,2,3,4,5], name=1)
    except Data.BadSpecifiedName: pass
    else: assert(False)

    # Try adding a column with multiple types
    try:
        b = a[:]
        b.add_column([1,2,3,4,5.0])
    except Data.BadValue: pass
    else: assert(False)

    # Try a mismatched-type slice assignment
    try:
        b = a[:]
        b[:,0] = [1.0, 1, 1, 1, 1]
    except Data.BadAssignment: pass
    else: assert(False)

    # Try a mismatched-type column assignment operation
    try:   a["0"] = list(a["0"])[:-1] + ["0"]
    except Data.BadAssignment: pass
    else: assert(False)

    # Try a too-short column assignment operation
    try:   a["0"] = list(map(str,a["0"]))[:-1]
    except Data.BadValue: pass
    else: assert(False)

    # Try a too-long column assignment operation
    try:   a["0"] = list(map(str,a["0"])) + [1]
    except Data.BadAssignment: pass
    else: assert(False)

    # Try a bad combination of names and types
    try:   Data(names=["a","b","c"], types=[str, str])
    except Data.NameTypeMismatch: pass
    else: assert(False)

    # Try a bad value in "names"
    try:   Data(names=[1,2,3], types=["a","b","c"])
    except Data.BadSpecifiedName: pass
    else: assert(False)

    # Try a bad value in "types"
    try:   Data(names=["a","b","c"], types=[1,2,3])
    except Data.BadSpecifiedType: pass
    else: assert(False)

    # Try bad value
    try:   a.append(["a","b","c"])
    except Data.BadValue: pass
    else: assert(False)

    # Try bad element
    try:   a.append([1,2])
    except Data.BadElement: pass
    else: assert(False)

    # Try non-existing column index
    try:   a["hello"]
    except Data.UnknownName: pass
    else: assert(False)

    # Try column indexing an uninitialized data
    try:   Data()["a"]
    except Data.Empty: pass
    else: assert(False)

    # Done testing
    print("passed.")

    # Convert the following into test cases for "Row" object.
    print()
    b = a[:]
    print(a)
    print(type(a.types), a.types)
    print(type(a.names), a.names)
    print(type(a[0]))
    print(a[0])
    print(a['0'])
    print(a[0][0])
    a[-1][1] = '1'
    print(a)
    a[-1][-1] = 1.
    print(a)
    a[-2][1] = None
    print(a)

# =============================================================
#      END OF Python 'Data' with named and typed columns     
# =============================================================


# Run test cases
if __name__ == "__main__":
    test_data()
