/* Program to perform 1D cellular automaton (CA) computations and to use 1D CA
   to solve the density classification problem.

  Skeleton program written by Artem Polyvyanyy, http://polyvyanyy.com/,
  September 2024, with the intention that it be modified by students
  to add functionality, as required by the assignment specification.
  All included code is (c) Copyright University of Melbourne, 2024.

  Authorship Declaration:

  (1) I certify that except for the code provided in the initial skeleton file,
  the program contained in this submission is completely my own individual
  work, except where explicitly noted by further comments that provide details
  otherwise. I understand that work that has been developed by another student,
  or by me in collaboration with other students, or by non-students as a result
  of request, solicitation, or payment, may not be submitted for assessment in
  this subject. I understand that submitting for assessment work developed by
  or in collaboration with other students or non-students constitutes Academic
  Misconduct, and may be penalized by mark deductions, or by other penalties
  determined via the University of Melbourne Academic Honesty Policy, as
  described at https://academicintegrity.unimelb.edu.au.

  (2) I also certify that I have not provided a copy of this work in either
  softcopy or hardcopy or any other form to any other student, and nor will I
  do so until after the marks are released. I understand that providing my work
  to other students, regardless of my intention or any undertakings made to me
  by that other student, is also Academic Misconduct.

  (3) I further understand that providing a copy of the assignment specification
  to any form of code authoring or assignment tutoring service, or drawing the
  attention of others to such services and code that may have been made
  available via such a service, may be regarded as Student General Misconduct
  (interfering with the teaching activities of the University and/or inciting
  others to commit Academic Misconduct). I understand that an allegation of
  Student General Misconduct may arise regardless of whether or not I personally
  make use of such solutions or sought benefit from such actions.

  Signed by: Dhruv Verma
  Dated:     29/09/2024
*/

#include <stdio.h>
#include <stdlib.h>
#include <assert.h>
#include <math.h>

/* #DEFINE'S -----------------------------------------------------------------*/
#define SDELIM "==STAGE %d============================\n"   // stage delimiter
#define MDELIM "-------------------------------------\n"    // delimiter of -'s
#define THEEND "==THE END============================\n"    // end message

#define NBRHD_SIZE 3    // maximum neighbourhood size
#define NBRHDS 8        // number of possible neighbourhoods

#define ONSTATEA 42     // signifies cell is on in ASCII
#define OFFSTATEA 46    // signifies cell is off in ASCII
#define ONSTATEN 1      // signifies cell is on in decimal
#define OFFSTATEN 0     // signifies cell is off in decimal
#define ONSTATEC '*'    // signifies cell is on as character
#define OFFSTATEC '.'   // signifies cell is on as character

#define NUMBER 0        // need to convert automaton to numbers
#define CHARACTER 1     // need to convert automaton to characters
#define RULE184 184     // update rule number for 184
#define RULE232 232     // update rule number for 232
#define COUNTCELL 2     // number of times specific cell needs to be counted
#define COUNT0 0        // access count cell array at index 0
#define COUNT1 1        // access count cell array at index 1

/* TYPE DEFINITIONS ----------------------------------------------------------*/
typedef struct {     // neighbourhood which corresponding state from update rule
    int neighbourhood[NBRHD_SIZE];     // neighbourhood of size three
    int state;                         // cell value of neighbourhood
} rule_t;

typedef struct {  // stores position and run for counting cell changes over time
    int cell_pos;      // cell position in automaton
    int start_run;     // starting run to count form
} count_t;

typedef struct {     // contains automaton sequence
    int *sequence;   // sequence
} state_t;

typedef struct {                     // elementary CA defined by
    int     size;                    // automaton size
    int     rule_num;                // update rule number
    int     runs;                    // number of evolutions
    rule_t  update_rule[NBRHDS];     // update rule with neighbourhood and state
    state_t *automaton;              // array of all automatons
    count_t cell_count[COUNTCELL];   // cells to count change over time
} CA_t;

/* USEFUL FUNCTIONS ----------------------------------------------------------*/
void stage0(CA_t *automata);
void stage1(CA_t *automata);
void stage2(CA_t *automata);
void read_input(CA_t *automata);
int runs_184(int size);
int runs_232(int size);
void convert_to_binary(int decimal, int binary_array[]);
void create_update_rule(rule_t update_rule[], int state_changes[]);
void print_update_rule(rule_t update_rule[]);
void convert_automaton(int automaton[], int size, int form);
void print_automaton(int automaton[], int size);
void update_state(CA_t *automata, int run, rule_t update_rule[]);
void update_cells_184(CA_t *automata);
void update_cells_232(CA_t *automata);
void cell_history(CA_t *automata, int end, int n);
void density_classification(CA_t *automata);
void free_memory(CA_t *automata);

/* WHERE IT ALL HAPPENS ------------------------------------------------------*/
/* main program provides traffic control */
int main(int argc, char *argv[]) {
    CA_t automata;
    read_input(&automata);
    stage0(&automata);
    stage1(&automata);
    stage2(&automata);
    free_memory(&automata);
    return EXIT_SUCCESS;        
}

/* USEFUL FUNCTIONS ----------------------------------------------------------*/
/* traffic control for stage0, creating and printing update rule and initial 
   automaton sequence */
void stage0(CA_t *automata) {
    // print size and rule from user input
    printf(SDELIM, 0);
    printf("SIZE: %d\nRULE: %d\n", automata->size, automata->rule_num);
    printf(MDELIM);

    // create and print update rule
    int state_changes[NBRHDS];
    convert_to_binary(automata->rule_num, state_changes);
    create_update_rule(automata->update_rule, state_changes);
    print_update_rule(automata->update_rule);
    printf(MDELIM);

    // print original automaton sequence
    printf("%4d: ", 0);
    print_automaton(automata->automaton[0].sequence, automata->size);
}

/* traffic control for stage1, updating and printing the automaton sequence 
   using user update rule for specified number of runs. specific cell position
   in automaton's states are counted as number of off and on states */
void stage1(CA_t *automata) {
    printf(SDELIM, 1);

    // print original automaton sequence
    printf("%4d: ", 0);
    print_automaton(automata->automaton[0].sequence, automata->size);

    // evolve automaton sequence for number for required runs and print each run
    for (int i = 1; i <= automata->runs; i++) {
        update_state(automata, i, automata->update_rule);
        printf("%4d: ", i);
        print_automaton(automata->automaton[i].sequence, automata->size);
    }
    printf(MDELIM);
    // determine number of on, off cells at specific cell position in automaton 
    cell_history(automata, automata->runs, COUNT0);
}

/* traffic control for stage2, solving automaton density problem by evolving
    automaton through update rule 184 and then rule 232 to find ratio of on and 
    off cells in original sequence. specific cell position in automaton's states
    are counted as number of off and on states */
void stage2(CA_t *automata) {
    // update automata using rule 184
    printf(SDELIM, 2);
    update_cells_184(automata);

    // update automata using rule 232
    printf(MDELIM);
    update_cells_232(automata);

    // determine number of on, off cells at specific cell position in automaton 
    printf(MDELIM);
    int total_runs =  automata->runs + runs_184(automata->size) + 
    runs_232(automata->size);
    cell_history(automata, total_runs, COUNT1);

    // print original automata 
    printf(MDELIM);
    printf("%4d: ", automata->runs);
    print_automaton(automata->automaton[automata->runs].sequence, 
                    automata->size);
    // density classification for original automata
    printf("AT T=%d: #ON/#CELLS ", automata->runs);
    density_classification(automata);
    printf(THEEND);
}

/* reads input form user into pointer to type CA_t from an input file */
void read_input(CA_t *automata) {
    // read automata size and update rule number
    scanf("%d\n%d\n", &(automata->size), &(automata->rule_num));
    
    // temporary character array to store starting automaton sequence as string
    char start_seq[automata->size + 1];
    scanf("%s\n", start_seq);
    
    // read automata runs, and cell position and starting run for cell counting
    scanf("%d\n%d,%d\n%d,%d\n", &(automata->runs), 
    &(automata->cell_count[0].cell_pos), &(automata->cell_count[0].start_run),
    &(automata->cell_count[1].cell_pos), &(automata->cell_count[1].start_run));

    // calculate the total runs required to store
    int total_runs = automata->runs + runs_184(automata->size) + 
    runs_232(automata->size) + 1;

    // assign memory of automaton to be large enough for total runs
    automata->automaton = (state_t *) malloc(sizeof(state_t) * total_runs);
    assert(automata->automaton != NULL);

    // assign memory of automaton sequence to be large enough for automata size
    for (int i = 0; i < total_runs; i++) {
        automata->automaton[i].sequence = (int *) malloc(sizeof(int) * 
        automata->size);
        assert(automata->automaton[i].sequence != NULL);
    }

    // convert initial into numbers 0, 1 and copy into first automaton sequence
    for (int i = 0; i < automata->size; i++) {
        if (start_seq[i] == ONSTATEC) {
            automata->automaton[0].sequence[i] = ONSTATEN;
        } else {
            automata->automaton[0].sequence[i] = OFFSTATEN;
        }
    }

}

/* calculates and returns number of runs required for rule 184 */
int runs_184(int size) {
    int runs = (size - 2) / 2;
    return runs;
}

/* calculates and returns number of runs required for rule 232 */
int runs_232(int size) {
    int runs = (size - 1) / 2;
    return runs;
}

/* takes update rule as integer and an integer array, converting the rule number
   into binary and stores the binary representation in the integer array */
void convert_to_binary(int decimal, int binary_array[]) {
    // array of exponents of binary
    int exponents[NBRHDS] = {7, 6, 5, 4, 3, 2, 1, 0};
    /* cycles through each exponent and if number is greater than binary 
       exponent, same position in binary array is assigned to be 1, otherwise 
       value is assigned to be 0  */
    for (int i = 0; i < NBRHDS; i++) {
        if (decimal >= pow(2, exponents[i])) {
            decimal -= pow(2, exponents[i]);
            binary_array[i] = ONSTATEN;
        } else {
            binary_array[i] = OFFSTATEN;
        }
    }
}

/* takes input array of type rule_t and state changes and creates an update rule
   based off all possible neighbourhood combinations and specifed cell changes 
   binary representation of update rule number */
void create_update_rule(rule_t update_rule[], int state_changes[]) {
    // array of all possible combinations of possible neighbourhoods
    int possible_nbrhds[NBRHDS][NBRHD_SIZE] = {{0,0,0}, {0,0,1}, {0,1,0}, 
    {0,1,1}, {1,0,0}, {1,0,1}, {1,1,0}, {1,1,1}};

    // each neighbourhood is assigned state value based off state_changes array
    for (int i = 0; i < NBRHDS; i++) {
        for (int j = 0; j < NBRHD_SIZE; j++) {
            update_rule[i].neighbourhood[j] = possible_nbrhds[i][j];
        }
        update_rule[i].state = state_changes[NBRHDS - 1 - i];
    }
}

/* prints update rule */
void print_update_rule(rule_t update_rule[]) {
    // printing all neighbourhoods in update rule
    for (int i = 0; i < NBRHDS; i++) {
        printf(" ");
        for (int j = 0; j < NBRHD_SIZE; j++) {
            printf("%d", update_rule[i].neighbourhood[j]);
        }
    }
    printf("\n");

    // printing all neighbourhood state values in update rule
    for (int i = 0; i < NBRHDS; i++) {
        printf("  %d ", update_rule[i].state);
    }
    printf("\n");
}

/* takes automaton sequence and converts it to either characters or integers  */
void convert_automaton(int sequence[], int size, int form) {
    /* if desired form is character, sequence is converted to integer value of 
       characters  '*' and '.' */
    if (form == CHARACTER) {
        for (int i = 0; i < size; i++) {
            if (sequence[i] == ONSTATEN) {
                sequence[i] = ONSTATEA;
            } else {
                sequence[i] = OFFSTATEA;
            }
        }
    } else { // desired form is number, sequence is converted to 0 and 1's
        for (int i = 0; i < size; i++) {
            if (sequence[i] == ONSTATEA) {
                sequence[i] = ONSTATEN;
            } else {
                sequence[i] = OFFSTATEN;
            }
        }
    }
}

/* prints automaton sequence */
void print_automaton(int sequence[], int size) {
    // converts sequence to character
    convert_automaton(sequence, size, CHARACTER);
    // cycle through sequence and print character
    for (int i = 0; i < size; i++) {
        printf("%c", sequence[i]);
    }
    printf("\n");
    // convert sequence back to numbers
    convert_automaton(sequence, size, NUMBER);
}

/* updates automaton sequence by based off update rule provided by comparing 
   left and right values of each cell position to all neighbourhoods and 
   assigning appropriate cell state within update rule */
void update_state(CA_t *automata, int run, rule_t update_rule[]) {
    int middle, left, right;
    // cycle through all cell positions of cell automata
    for (int i = 0; i < automata->size; i++) {
        middle = automata->automaton[run - 1].sequence[i];

        // assign left and right cell positions wrapping around automata
        if (i == 0) {
            left = automata->automaton[run - 1].sequence[automata->size - 1];
            right = automata->automaton[run - 1].sequence[i + 1];
        } else if (i == automata->size - 1) {
            left = automata->automaton[run - 1].sequence[i - 1];
            right = automata->automaton[run - 1].sequence[0];
        } else {
            left = automata->automaton[run - 1].sequence[i - 1];
            right = automata->automaton[run - 1].sequence[i + 1];
        }

        /* cycle through each update rule and check if set of cells match and 
           change state of middle cell once match is found */
        for (int j = 0; j < NBRHDS; j++) {
            if (left == update_rule[j].neighbourhood[0] && 
                middle == update_rule[j].neighbourhood[1] &&
                right == update_rule[j].neighbourhood[2]) {
                    automata->automaton[run].sequence[i] = 
                    update_rule[j].state;
                    break;
                }
        }
    }
}

/* creates update rule 184 and updates the automaton sequence based off rule 184
   based of required runs for rule 184 */
void update_cells_184(CA_t *automata) {
    // calculate number of runs required to process rule 184
    int run_184 = runs_184(automata->size);
    printf("RULE: %d; STEPS: %d.\n", RULE184, run_184);

    // create update rule 184
    int state_changes[NBRHDS];
    rule_t rule_184[NBRHDS];
    convert_to_binary(RULE184, state_changes);
    create_update_rule(rule_184, state_changes);
    printf(MDELIM);

    // print inital automaton sequence
    printf("%4d: ", automata->runs);
    print_automaton(automata->automaton[automata->runs].sequence, 
                    automata->size);

    // update and print sequence based off number of runs calculated
    for (int i = 1; i <= run_184; i++) {
        update_state(automata, automata->runs + i, rule_184);
        printf("%4d: ", i + automata->runs);
        print_automaton(automata->automaton[automata->runs + i].sequence, 
                        automata->size);
    }
}

/* creates update rule 232 and updates the automaton sequence based off rule 232
   based of required runs for rule 232 */
void update_cells_232(CA_t *automata) {
    // calculate number of runs required to process rule 232
    int run_232 = runs_232(automata->size);
    int run_184 = runs_184(automata->size);
    printf("RULE: %d; STEPS: %d.\n", RULE232, run_232);
    
    // create update rule 232
    int state_changes[NBRHDS];
    rule_t rule_232[NBRHDS];
    convert_to_binary(RULE232, state_changes);
    create_update_rule(rule_232, state_changes);
    printf(MDELIM);

    // print initial automata
    printf("%4d: ", automata->runs + run_184);
    print_automaton(automata->automaton[automata->runs + run_184].sequence, 
                    automata->size);

    // update and print sequence based off number of runs calculated
    for (int i = 1; i <= run_232; i++) {
        update_state(automata, automata->runs + run_184 + i, rule_232);
        printf("%4d: ", i + automata->runs + run_184);
        print_automaton(automata->automaton[automata->runs + 
                        run_184 + i].sequence, automata->size);
    }
}

/* counts the number of on vs off states of cell in specific cell position in 
   automata and prints number of off and on cells at that position starting from
   specific run */
void cell_history(CA_t *automata, int end, int n) {
    // identify starting position to count cells from
    int pos = automata->cell_count[n].cell_pos;
    int off = 0, on = 0;

    /* cycle through sequence checking whether state is on and off, incrementing
       respective counters */
    for (int i = automata->cell_count[n].start_run; i <= end; i++) {
        int state = automata->automaton[i].sequence[pos];

        if (state == ONSTATEN) {
            on = on + 1;
        } else {
            off = off + 1;
        }
    }
    // print number of on and off cells
    printf("#ON=%d #OFF=%d CELL#%d START@%d\n", on, off, 
    automata->cell_count[n].cell_pos, automata->cell_count[n].start_run);
}

/* solves density classification program by identifying whether last automaton 
   contains all off, all on or alternating on and off cells and prints whether
   original automaton has more off or on states in its sequence */
void density_classification(CA_t *automata) {
    // calculate total runs of automaton evolutions
    int total_runs =  automata->runs + runs_184(automata->size) + 
    runs_232(automata->size);
    int on = 0, off = 0;

    // count number of on and off cells in final automaton sequence
    for (int i = 0; i < automata->size; i++) {
        if (automata->automaton[total_runs].sequence[i] == ONSTATEN) {
            on++;
        } else {
            off++;
        }
    }

    // classify density of automaton based on ratio of on vs off cells
    if (on > off) {
        printf("> 1/2\n");
    } else if (on < off) {
        printf("< 1/2\n");
    } else {
        printf("= 1/2\n");
    }
}

/* frees memory manually allocated to automata */
void free_memory(CA_t *automata) {
    // calculates total size of automaton array
    int total_runs = automata->runs + runs_184(automata->size) + 
    runs_232(automata->size) + 1;
    
    // cycles through each  sequence and frees memory, assigning null pointer
    for (int i = 0; i < total_runs; i++) {
        free(automata->automaton[i].sequence);
        automata->automaton[i].sequence = NULL;
    }
    // frees automaton array, assigning null pointer
    free(automata->automaton);
    automata->automaton = NULL;
}

/* THE END -------------------------------------------------------------------*/
// algorithms are fun