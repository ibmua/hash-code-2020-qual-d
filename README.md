# My Algorithm to solve test D
⚠ SPOILER ALERT

Got max-25 books on D. It was quite tough. 
Solution is this

### REPEAT
➊ Put libs that can yield +books back on, prioritizing those libs that add most books.

➋ Remove some libs - first those removing which doesn't remove any books, then those that are only required for 1 book, then 2.

_I first tried removal in a greedy way with heuristic that prioritizes removal of libs that will disable the least number of other libs from being removed. Don't bother.._

__Add randomization of removal and addition sequences and repeat over 10000-100000 iterations.__ For example, when adding books, if two books can yield same amount of +unique_books we can randomly choose which one of them to add first. When removing books you can actually ignore and don't bother much with heuristics that I _wrote above_. __Randomization and repeating is key.__ As far as I tested, if not repeating, the best heuristic I made, which took like an hour+ to figure and test left me with a little over 15500 libs and all the books. While keeping heuristics at a minimum and just injecting random and repeating yields results like 15025 while keeping all of the books. Getting to 15500 will take just a few add-subtract cycles. So randomizing and repeating is key here. Heuristics - don't bother.

I repeat that for many many iterations and slowly the worst libs sort of fall off while the best ones crystallize by becoming too key to be peeled off. They are key to more books than those which we peel off.
So for me the amount of libs that provide any books went down from >15500 to 15025.
Then I just peel off the libs that provide me with only a single unique book and put them at the end of the output sequence, while the rest is put at the beginning.

I think this is a great example of __antifrigility__.
