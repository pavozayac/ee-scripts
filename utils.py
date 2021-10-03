def chunks(lst: list, n: int):
       for i in range(0, len(lst), n):
              yield lst[i:i + n]