import pytest
from functions import *


def test_reverse_string():
    assert reverse_string("hello") == "olleh"
    assert reverse_string("Python") == "nohtyP"
    assert reverse_string("") == ""


@pytest.mark.parametrize(
    "num, expected", [(2, True), (17, True), (15, False), (1, False), (0, False)]
)
def test_is_prime(num, expected):
    assert is_prime(num) == expected


@pytest.fixture
def default_user():
    return User(username="test_user", email="test_user@example.com")


def test_user_username(default_user):
    assert default_user.username == "test_user"


def test_user_email(default_user):
    assert default_user.email == "test_user@example.com"


def test_divide():
    assert divide(10, 2) == 5
    with pytest.raises(ZeroDivisionError):
        divide(1, 0)
