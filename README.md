Dear Fellow Participant,

It seems that you have been shared a very curious repository. 

Many years of coding experience should tell you that this repository, named after the famous game Othello (or, if you peasants prefer, Reversi, but I will use Othello as it sounds classy), is related to Othello. Unlike the tests at UC Berkeley, there are no tricks. So yes, this repository is about Othello.

And it seems that you have been chosen as a collaborator, or a spectator (shade).

Currently, the repository contains a folder called "FinalSolver." It is called FinalSolver not because its code is finalized, but instead it is final in its algorithmic design, and it is also final because it works. This means I have not made that many mistakes in my code. Round of applause for me. 

Notice that the "MPISolver" folder has one file, and the file is empty. That is in fact the very next task. The goal is to implement an Othello solver with UC Berkeley's supercomputer Savio. Unfortunately, that requires a lot more work, and therefore it is still an empty folder, desperately waiting for meaning to be added into its emptiness.

After all this nonsense, let us get into the serious business:

First, please search up the rules of Othello if you do not know them. The point of this repository is to solve the game, not to play it, so we won't be explaining the rules here. 

The overarching theme of this solver is an attempt to "parallelize" solving. There are obviously many attempts in the past, but a recent attempt by an alumni, Justin Yokota, used the idea of "shards" to solve Connect 4. The idea was that, if we created "shards" of the game tree of Connect 4 by using a column (in his paper, he uses the leftmost column), then we can effectively create a game tree with condensed nodes and little to no edges between nodes. After all, the "lower" pieces of a particular column can never change. In this way, he was able to create sections of the game tree that can be solved in parallel. We use a similar and more generalized (in some ways) approach. Feel free to read his paper at this link: https://www2.eecs.berkeley.edu/Pubs/TechRpts/2022/EECS-2022-219.pdf

If you read Justin's paper, one proposal to "shard" Othello by using the four corners of the game board, as those pieces will never change color. In that way, we are able to break these boards apart into shards. He extends this to pieces adjacent to the corners (roughly speaking) as they will also not change in color. This can create finer shards. However, he also notes that this type of hashing algorithm will not be monotonic and will be difficult to deal with.

Therefore, we approach this differently. Notice that in order to parallelize solving a board game, we want to split the game tree into similar sized nodes, and we want it to have, ideally, no edges between these nodes. We generalize this idea to, what I call, a "blob."

A "blob" is an abstract idea similar to a shard. However, it is different as instead of taking "values" of pieces (white or black in Othello, red or yellow in Connect 4), it takes some sort of "metadata" of the board to split up the game tree. This allows a higher degree of generalization. The reason why we call it a "blob" is because of our initial use case. In Othello, we define a "blob" to be a board state where pieces are treated all the same. In other words, for a given board, we can define a "blob" by placing 0s in positions where there is no piece, and 1s in positions where there is a piece. This allows us to create a hash value for our "blob" by simply taking the face value of this number. 

We then use the idea of page tables. Notice that each "blob" is similar to a "page" due to the fact that the range of hashes is roughly fixed. If we represent each piece as 0 for black and 1 for white, then we have that the maximum hash is the hash of all 1s. In other words, if all the pieces are white. In Othello, for example, this value has a maximum value of 2**64 - 1 for the maximum blob, a very beautiful number. If two "blobs" are of the same size, we have that they have the same maximum values. This means that, roughly speaking, most of these "blobs" will be about the same size for a given tier, tier being defined as the number of pieces on the board (which is different from remoteness, as in Othello you may be forced to skip turns). 

Combining these ideas, we are able to traverse the game tree while distributing work fairly equally to workers. We start with an initial position, and during discovery phase, we place the boards (after accounting for symmetries) with the same "blob" together in a file. For each specific board, we define black pieces as 0, white pieces as 1, and we go through the pieces in order (right to left, bottom to top, or in any way you prefer. We do it in this way since we read numbers from left to right, top to bottom), ignoring empty spaces. This gives us a "blob" specific hash, and we store the FOUND macro in this position. Then, we continue this process until we reach the end of the game tree, keeping track of the order of "blobs" discovered (this will be important later).

Then, reversing the order of discovered "blobs," we start distributing work of same-tiered blobs (remember, we define tier as the number of pieces now, not remoteness), and we calculate remoteness and primitive values starting from the bottom. After discovering a tier, we move to the next tier, which starts at whatever is next and is contiguous as we store it during the discovery phase in order. This continues until we return to the initial position. 

After this, we are done.

This method also allows for fast retrieval. For any given position (after symmetries), we first determine the "blob," then determine the "blob" specific hash. Then, we open the "blob" file and go to that hash position to retrieve the remoteness and primitive value of that position. 

We also implement efficient compression. This is directly stolen from Justin Yokota's paper, where we abuse the efficient compression of repetitive information with gzip to decrease the zipped database file. 

All in all, our method decreases the solving time (compared to the GamesCrafter's general solver) by quite a bit, as solving it in this method uses about 1 second, whereas the general solver uses numerous minutes. Furthermore, the compression of the database is significant. Comparing to GamesCrafter's general solver, we decrease the database significantly even before compression, with a total size of 842 kb. After compression with 7zip (as this compresses across files, which is what we want for fast retrieval), we decrease the database to a mere 13.6 kb, which is about 1.6% of our decompressed database. The number of discovered positions is 10513 canonical positions, so about 11 bits per position after compression. Honestly, I do not know how to improve this (or if it is even possible), so any suggestions would be great. 

All of the searching here assumes that we maximize symmetries. We have 8 symmetries for each board, and we always find the one with the smallest "blob" hash, then smallest "blob" specific hash. There are also moves that are "skip" moves because one player cannot make a move. We only store the position where the player could make a move, the "skip" moves can be extrapolated from those moves. If you ask what happens if both players couldn't make a move. Well, that is a primitive position, and you would have known that if you read the rules of Othello. So if you didn't, please go read the rules again, then reread this README to pay for your sin of laziness.


To compile the one core code, please cd into FinalSolver, then run:
gcc -o solver solver.c othello4.c file_lst.c num_lst.c game.c

Then, you can run:
./solver

This will generate Data4x4 (assuming you pulled it straight out of the repository). 

This is the end of the naive description. I hope you enjoy working on or viewing this project. Please contact Abraham Hsu at hpghsu@berkeley.edu or hpghsu@gmail.com for any inquiries. 