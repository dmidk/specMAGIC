# Coding standard

Below is a list of good-practice guidelines for use when contributing to specMAGIC. 

## Copy constructors 

Any struct or class in which you expect to store large quantities of data - i.e. anything more than a tiny struct - should have its copy constructor and copy assignment operator explicitly deleted (Singleton implementation). This is to avoid accidental copies of large amounts of data, which kills performance. See the `Matrix` struct, which deletes the copy constructor and assignment operator and explicitly defines the move constructor and assignment operator, for good practice. 

## The `<vector>` library

The basic C-style array functionality is not particularly convenient or flexible. It is not possible to resize on-the-fly, for example. The most popular solution to this problem is the `<vector>` library. 

`vector` provides significant flexibility and behaves much as the `numpy` arrays common to `python` codes. However, the shared memory constructs used in the parallelised portion of the code require careful handling. It is very easy to introduce race conditions and nondeterministic behaviour if the array holding some relevant data may be being dynamically resized without the user being aware.

For this reason, **do not** use the `<vector>` class in any performance-intensive part of the program. No `vector` object should appear inside any `OpenMP` parallel region. It is perfectly fine to use it in e.g. read-in or write-out routines, where the `vector` objects are temporary and will not be carried into any parallel loop. 

There are many instances where `vector` works well with OpenMP; this rule is overcautious. For anything performance-critical, prefer the method used for the `Matrix` struct, which contains a flat C-style array with accessor functions built on top for convenience. This is fast, cheap access which throws an immediate error if the program attempts to access illegal memory. 

## Variable typing 

Custom types used in the program are defined in `types.hpp`. Furthermore, avoid the use of `size_t` (as its actual size is hardware-dependent) and prefer `int`, `unsigned int`, `long` or similar as appropriate. 

Unfortunately `size_t` is required for some `netcdf` operations, as this library requires it. This is unavoidable. Do not use it elsewhere.

## Others

In general, specMAGIC uses the `i++` rather than `++i` style for loops. 

specMAGIC uses `struct` object rather than `class`. Use a `struct` rather than a `class`, unless you need members to be private by default, for consistency.

Be wary of including large libraries in the global namespace, as this can affect performance. The `std` library, for example, is huge. Avoid `using namespace std` and similar wherever possible. 

This list is subject to ongoing updates.