# class used as function "decorator":
class Memoized:
  def __init__(self, function):
    self._function = function
    self._cache = {}
  def __call__(self, *args):
    if args not in self._cache:
      # not in the cache: call the function and store the result in
      # the cache
      self._cache[args] = self._function(*args)
    # the result must now be in the cache:
    return self._cache[args]

