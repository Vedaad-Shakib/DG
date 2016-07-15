import re

test_template = open("xf_Math.test.c", "r").read()
test_write = open("xf_Math.test.transformed.c", "w")

# find and replace function names
test_template = re.sub(r'TEST_([a-zA-Z1-9_]+)\(\)', r'TEST(\1, All)', test_template)
# find and eliminate return
test_template = test_template.replace("return xf_OK;", "")
# find and replace assert functions
test_template = re.sub(r'xf_AssertRealVectorWithin\(([a-zA-Z0-9_\-\.*]+), ([a-zA-Z0-9_\-\.*]+), ([a-zA-Z0-9_\-\.*]+), ([a-zA-Z0-9_\-\.*]+)\);', r'EXPECT_TRUE(AssertRealVectorWithin(\1, \2, \3, \4));', test_template);

test_template = re.sub(r'xf_AssertEqual\(([a-zA-Z0-9_\-\.*]+), ([a-zA-Z0-9_\-\.*]+)\);', r'EXPECT_EQ(\1, \2);', test_template);

test_template = re.sub(r'xf_AssertWithin\(([a-zA-Z0-9_\-\.*]+), ([a-zA-Z0-9_\-\.*]+), ([a-zA-Z0-9_\-\.*]+)\);', r'EXPECT_DOUBLE_EQ(\1, \2);', test_template);

test_template = re.sub(r'xf_AssertIntVectorEqual\(([a-zA-Z0-9_\-\.*]+), ([a-zA-Z0-9_\-\.*]+), ([a-zA-Z0-9_\-\.*]+)\);', r'EXPECT_TRUE(AssertIntVectorEqual(\1, \2, \3));', test_template);

# find and eliminate comments and includes
test_template = re.sub(r'/\*.*?\*/\n*', '', test_template)
#test_template = re.sub(r'#include.*\n', '', test_template)

test_write.write(str(test_template))
