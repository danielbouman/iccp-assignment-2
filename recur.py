
def factorial(n,m):
    if n == 1:
        print(m)
        return 1
    else:
        print('n!=1 '+str(m))
        return n * factorial(n-1,m)

result = factorial(4,5)
print(result)