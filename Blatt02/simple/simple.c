/*
** simple error demonstration to demonstrate power of valgrind
** Julian M. Kunkel - 17.04.2008
*/

#include <stdio.h>
#include <stdlib.h>

  
/*
 * allocated memory space at the adress before writing the array, thereby ensuring the local variable is not overwritten.
 */
int *mistake1(void) {
  int *buf = malloc(sizeof(int) * 6);
  buf[0] = 1;
  buf[1] = 1;
  buf[2] = 2;
  buf[3] = 3;
  buf[4] = 4;
  buf[5] = 5;
  return buf;
}


/*
 * changed "sizeof(char) * 4" to "sizeof(int) * 4" as chars are not big enough. buf[2] would thereby be accessing out of bounds.
 * changed "buf[2] = 2" to "buff[1] =2".
 */
int *mistake2(void) {
  int *buf = malloc(sizeof(int) * 4);
  buf[1] = 2;
  return buf;
}

int *mistake3(void) {
  /* In dieser Funktion darf kein Speicher direkt allokiert werden. */
  int *buf = (int *) mistake2();
  buf[0] = 3;
  return buf;
} 

//mistake 2

/*
 * changed index writing to buf, as main() calls the first element.
 * removed "free()"-call before accessing the array
 */
int *mistake4(void) {
  int *buf = malloc(sizeof(int) * 4);
  buf[0] = 4;
  //free(buf);
  return buf;
}

/**
 * no need to free p[3], as no malloc() happens. 
 * 
*/

int main(void) {
  /* Modifizieren Sie diese Zeile nicht! */
  int *p[4] = {&mistake1()[1], &mistake2()[1], mistake3(), mistake4()};

  printf("1: %d\n", *p[0]);
  printf("2: %d\n", *p[1]);
  printf("3: %d\n", *p[2]);
  printf("4: %d\n", *p[3]);

  /* mhh muss hier noch etwas gefreed werden? */
  /* Fügen sie hier die korrekten aufrufe von free() ein */
  free(p[3]); /* welcher Pointer war das doch gleich?, TODO: Fixme... ;-) */
  free 
  return 0;
}
