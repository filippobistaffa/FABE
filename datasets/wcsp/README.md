Weighted Constraint Satisfaction Problems (WCSP)
===================
WCSP is a text format composed of a list of numerical and string terms separated by spaces. Instead of using names for making reference to variables, variable indexes are employed. The same for domain values. All indexes start at zero.

Cost functions can be defined in intention (see below) or in extension, by their list of tuples. A default cost value is defined per function in order to reduce the size of the list. Only tuples with a different cost value should be given. All the cost values must be positive. The arity of a cost function in extension may be equal to zero. In this case, there is no tuples and the default cost value is added to the cost of any solution. This can be used to represent a global lower bound of the problem.

Definition of Problem Name, Problem Dimensions, and Domain sizes
----------
The first line contains the following information:

    <Problem name> 
    <Number of variables (N)>
    <Maximum domain size>
    <Number of cost functions>
    <Initial global upper bound of the problem (UB)>
 
 The second line contains the sizes of the domains of all variables:

    <Domain size of variable with index 0>
    ...
    <Domain size of variable with index N-1>

Definition of Cost Functions
----------
Definition of a cost function in extension:

    <Arity of the cost function>
    <Index of the first variable in the scope of the cost function>
    ...
    <Index of the last variable in the scope of the cost function>
    <Default cost value>
    <Number of tuples with a cost different than the default cost>

followed by for every tuple with a cost different than the default cost:

    <Index of the value assigned to the first variable in the scope>
    ...
    <Index of the value assigned to the last variable in the scope>
    <Cost of the tuple>

Definition of a cost function in intension by giving its keyword name and `K` parameters (and replacing the default cost value by -1):

    <Arity of the cost function>
    <Index of the first variable in the scope of the cost function>
    ...
    <Index of the last variable in the scope of the cost function>
    -1
    <keyword>
    <parameter1>
    ...
    <parameterK>
