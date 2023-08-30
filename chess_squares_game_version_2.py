# -*- coding: utf-8 -*-

class Square:
    
    '''
    A class that stores the basic info about a square from row and file number, 
    including color and coordinates
    '''
        
    def __init__(self, file_num, row):
        
        files = ['a', 'b', 'c', 'd', 'e', 'f', 'g', 'h']
        
        self.file_num = file_num
        self.row = row
        self.file = files[self.file_num - 1]
        self.square_coords = self.file + str(self.row)
    
        if (self.file_num % 2 == 0) & (self.row % 2 == 0):
            self.color = 'dark'
        elif (self.file_num % 2 != 0) & (self.row % 2 == 0):
            self.color = 'light'
        elif (self.file_num % 2 == 0) & (self.row % 2 != 0):
            self.color = 'light'
        elif (self.file_num % 2 != 0) & (self.row % 2 != 0):
            self.color = 'dark'
        

def getSquare():
    
    '''
    Outputs a random square on the chess board
    '''
    
    import random
    square = Square(random.randrange(1, 9), random.randrange(1, 9))
    
    return square

    
def colorQuiz():
    
    '''
    Loops through the game however many times are indicated by the player, 
    then returns the score at the end. 
    '''
    
    print('\nType exit at any time to quit\n')
    iterations = int(input('How many times would you like to play?\n'))
    print('\n')
    
    correct_answers = 0
        
    for i in range(iterations):
        
        square = getSquare()
        
        answer = input('Is ' + square.square_coords + ' a light square or a dark square?\n')
        
        if answer.lower() == 'exit':
            break
        elif answer.lower() == square.color or answer.lower() == square.color[0]: #allow for 'l' and 'd'
            print('\nCorrect!\n')
            correct_answers += 1
        else:
            print('\nIncorrect\n')
            
    print('\nFinal score: ' + str(correct_answers) + '/' + str(iterations))
    
    
def getLegalKnightMoves(starting_square):
    
    '''
    Takes a Square object as input and outputs all the legal moves a knight 
    can make as a dictionary with the name of the square as the key and 
    the Square object as the value
    '''
    
    movement_patterns = [[2, 1], [1, 2], [-2, 1], [-1, 2], [2, -1], [1, -2], [-2, -1], [-1, -2]]
    
    legal_moves = {}
    #legal_moves_coords = []
    
    for movement in movement_patterns: #go through the movement patterns from the starting point and record legal moves in legal_moves
        try:
            new_square = Square(starting_square.file_num + movement[0], starting_square.row + movement[1])
            if new_square.file_num > 0 and new_square.file_num < 9 and new_square.row > 0 and new_square.row < 9:
                legal_moves[new_square.square_coords] = new_square
        except IndexError:
            pass
    
    return legal_moves

    
def knightWalk():
    
    '''
    Allows you to walk a knight around the board, either from a randomly 
    generated square or from one of the four starting knight squares
    '''
    
    import random
    
    print('\nType exit at any time to quit\n')
    
    starting_square_preference = input('Would you like to start on a random square (1) or a normal knight starting square (2)?\n')
    if starting_square_preference == '1':
        current_square = getSquare()
    elif starting_square_preference == '2':
        knight_starting_squares = [Square(2, 1), Square(7, 1), Square(2, 8), Square(7, 8)]
        current_square = random.choice(knight_starting_squares)
    else:
        print('Idk what that means. I\'ll just give you a random square.')
        current_square = getSquare()
    iterations = int(input('How many steps would you like to take?\n'))

    steps = 0
    
    while steps < iterations:
        legal_moves = getLegalKnightMoves(current_square)
        
        print(f'\nThe knight is currently on {current_square.square_coords}')
        user_move = input('Enter your knight move\n')   
        
        if user_move.lower() == 'exit':
            break
        elif user_move in legal_moves.keys():
            current_square = legal_moves[user_move]
            print('\nThis is a legal move!')
            steps += 1
        else:
            print('\nGO FUCK YOURSELF DUMBASS.')

            
def chooseGame():
    
    player_choice = input('\nWhich game would you like to play? Press 1 or 2\n1. Color quiz\n2. Knight Walk\n')
    
    if player_choice.lower() == '1':
        colorQuiz()
    elif player_choice.lower() == '2':
        knightWalk()
    #else:
        #print('\nWhich game would you like to play? Enter 1 or 2\n1. Color quiz\n2. Knight Walk\n')


if __name__ == "__main__":
    
    chooseGame()
    
    play_again = input('\nWould you like to play again?\n')
    print("Your input is " + play_again.lower())
    while True:
        if play_again.lower() == ('yes' or 'y'):
            print('here is a fun spot')
            chooseGame()
            print('how about here')
        else:
            print('we are actually here')
            break
    
    print('We are here')
    
    
#%%
n = 0 
while True:
    
    if n < 4:
        n += 1
        print(n)    
    else: 
        break
    
