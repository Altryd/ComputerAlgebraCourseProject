## Мини-курсовая работа по предмету "Компьютерная алгебра"
### Алгоритм Полларда ECDLP

main.py - реализация всего что можно, включая ро-алгоритм полларда и нахождения порядка точки, на чистом питоне (почти) без sagemath/sympy/т.д.

unittests_pollard.py - тесты для main.py

sage_math.py - реализация ро-алгоритма Полларда и нахождение порядка точки на кривой


Проблемы:
- очень долго работает на больших числах

Для решения проблемы нужно посмотреть комментарии в sage_math.py для функций get_order_of_point(...) и func_(...).
Еще возможно лучше генерировать в func_(...) s штук Mi точек и их добавлять в зависимости от того, в какую область попал x.