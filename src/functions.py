def reverse_string(s):
    return s[::-1]

def is_prime(n):
    if n <= 1:
        return False
    for i in range(2, n):
        if n % i == 0:
            return False
    return True

class User:
  def __init__(self, username, email):
    self.username = username
    self.email = email

def divide(x, y):
    if y == 0:
        raise ZeroDivisionError("division by zero")
    return x / y