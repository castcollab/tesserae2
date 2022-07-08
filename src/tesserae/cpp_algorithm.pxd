# distutils: language=c++

cdef extern from "<algorithm>" namespace "std" nogil:
    # http://en.cppreference.com/w/cpp/algorithm/transform
    OutputIter transform[InputIter, OutputIter, UnaryOperation](InputIter, InputIter, OutputIter, UnaryOperation)

    # https://en.cppreference.com/w/cpp/algorithm/max_element
    ForwardIt max_element[ForwardIt](ForwardIt, ForwardIt)

    InputIt find[InputIt, T](InputIt, InputIt, T)


cdef extern from "<iterator>" namespace "std" nogil:
    long distance[InputIter](InputIter, InputIter)
