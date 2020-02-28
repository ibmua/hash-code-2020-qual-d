# My Algorithm to solve test D
⚠ SPOILER ALERT

Got max-47 books on D. It was quite tough. 
Solution is this

🔘 first of all add all libs

### REPEAT
➊ Remove some libs - first those removing which doesn't remove any books, then those that are only required for one book

➋ Put libs that can yield +books back on, prioritizing those libs that add most books.

__I do removal in a wise greedy way with heuristic that prioritizes removal of libs that will disable least number of other libs from being removed.__

__Some randomization is also added to greediness.__ For example when adding books, if two books can yield same amount of +unique_books we can randomly choose which one will be added first. __Randomization and repeating is key.__ As far as I tested, if not repeating, the best heuristic (without random) left me with a little over 15500 libs and all the books. While keeping heuristics at a minimum and just injecting random and repeating yields results like 15047 eventually while getting to 15500 will take just a few add-subtract cycles. So randomizing and repeating is key here. Heuristics not so much.

I repeat that for many many iterations and slowly the worst libs sort of fall off while the best ones crystallize by becoming too key to be peeled off. They are key to more books than those which we peel off.
So for me the amount of libs that provide any books went down from >15500 to 15047.
Then I just peel off the libs that provide me with only a single unique book and put them at the end of the output sequence, while the rest is put at the beginning.

*D – tough choices
5,105,945 points*
