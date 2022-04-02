 /** \file barnes.c
       \author Federico Campo
     Si dichiara che il contenuto di questo file e' in ogni sua parte opera
     originale dell'autore.  */

#include<stdio.h>
#include<string.h>
#include<stdlib.h>
#include<math.h>
double _s;

typedef struct nodo {
  /** coordinate e massa del corpo (foglia) 
      o del centro di massa (nodo interno) */
  double x, y; 
  double massa;

  // puntatori ai nodi figli 
  struct nodo * NW, *NE, *SE, *SW;

} nodo_t;

int stringa_a_corpo (const char* s, double* x, double* y, double* m){

  sscanf(s, "%lf %lf %lf", x, y, m);

  if(x != NULL && y != NULL && m !=NULL){
    return 0;}
  else return -1;
}


/* parametri della funzione: m, x, y, e root2 sono equivalenti a quelli della funzione inserisci, x_c e y_c 
sono le coordinate del centro del quadrante che di volta in volta vado a considerare e k è un parametro 
che viene dimezzato ad ogni chiamata ricorsiva della funzione che viene utilizzato per modificare il 
centro del quadrante da considerare
*/

nodo_t * inserisci_ricors(double m, double x, double y, nodo_t * root2, double x_c, double y_c, double k){

//dimezzo la variabile k che servirà per modificare di volta in volta le coordinate del centro
  k=k/2;

//creo i 4 figli e li pongo tutti con massa 0 e coordinate (0,0)
  if(root2->NE==NULL){
    root2->NE=malloc(sizeof(nodo_t)); if (root2->NE == NULL){perror("malloc"); exit(EXIT_FAILURE);}
    root2->NE->massa=0;
    root2->NE->x=0;
    root2->NE->y=0;}
  if(root2->SE==NULL){
	root2->SE=malloc(sizeof(nodo_t)); if (root2->SE == NULL){perror("malloc"); exit(EXIT_FAILURE);}
    root2->SE->massa=0;
    root2->SE->x=0;
    root2->SE->y=0;}
  if(root2->NW==NULL){
	root2->NW=malloc(sizeof(nodo_t)); if (root2->NW == NULL){perror("malloc"); exit(EXIT_FAILURE);}
    root2->NW->massa=0;
    root2->NW->x=0;
    root2->NW->y=0;}
  if(root2->SW==NULL){
    root2->SW=malloc(sizeof(nodo_t)); if (root2->SW == NULL){perror("malloc"); exit(EXIT_FAILURE);}
    root2->SW->massa=0;
    root2->SW->x=0;
    root2->SW->y=0;

//sposto il corpo che era prima nel nodo in uno dei suoi 4 figli, in base alla posizione in cui si trova rispetto alle coordinate del centro
     if(root2->x > x_c){
	
    if(root2->y >= y_c){
	root2->NE->massa=root2->massa;
      root2->NE->x=root2->x;
      root2->NE->y=root2->y;}
    else{//if y<y_c
	root2->SE->massa=root2->massa;
      root2->SE->x=root2->x;
      root2->SE->y=root2->y;}

  }

  else {	// if root2->x < x_c

    if(root2->y >= y_c)
     {root2->NW->massa=root2->massa;
      root2->NW->x=root2->x;
      root2->NW->y=root2->y;}
    else{ //if y<y_c
      root2->SW->massa=root2->massa;
      root2->SW->x=root2->x;
      root2->SW->y=root2->y;}
      }
  }





//cerco il punto in cui inserire il corpo, in base alle sue coordinate, confrontate con quelle del centro del quadrante o sottoquadrante che sto considerando 
  if (x > x_c){    
    if(y >= y_c){
  //Inserisco l'altro corpo se nel nodo non c'è nessun altro corpo (e quindi ha massa=0)
      if(root2->NE->massa==0){
	root2->NE->massa=m;
	root2->NE->x=x;
	root2->NE->y=y;
      }
//altrimenti richiamo la funzione ricorsivamente e modifico di volta in volta le coordinate del centro del quadrante o sottoquadrante
      else {root2->NE = inserisci_ricors(m, x, y, root2->NE, x_c + k  , y_c + k, k );}
    }
    
    else{  // y < y_c
      if(root2->SE->massa==0){

	root2->SE->massa=m;
	root2->SE->x=x;
	root2->SE->y=y;

      }
      else {root2->SE = inserisci_ricors(m, x, y, root2->SE , x_c + k , y_c - k, k );}
    }
  }

  
  else{    //se x < x_c
    if(y>= y_c){
      if(root2->NW->massa==0){
	root2->NW->massa=m;
	root2->NW->x=x;
	root2->NW->y=y;
      }
      else{root2->NW =inserisci_ricors(m, x, y, root2->NW, x_c - k , y_c + k, k );}
    }
    

    else{ //if y<y_c
      if(root2->SW->massa==0){
	root2->SW->massa=m;
	root2->SW->x=x;
	root2->SW->y=y;
      }
      else {root2->SW = inserisci_ricors(m, x, y, root2->SW, x_c - k, y_c - k, k);}
    }
  }

  //Aggiorno il nodo con le medie delle posizioni e lo faccio diventare centro di massa.

  root2->massa = (root2->NE->massa + root2->SE->massa + root2->NW->massa + root2->SW->massa);
  root2->x = ((root2->NE->massa)*(root2->NE->x) + (root2->SE->massa)*(root2->SE->x) + (root2->NW->massa)*(root2->NW->x) + (root2->SW->massa)*(root2->SW->x)) / (root2->massa);
  root2->y = ((root2->NE->massa)*(root2->NE->y) + (root2->SE->massa)*(root2->SE->y) + (root2->NW->massa)*(root2->NW->y) + (root2->SW->massa)*(root2->SW->y)) / (root2->massa);

  return root2;
}











nodo_t * inserisci(double m, double x, double y, nodo_t *root){


  //caso base, creo la radice principale e sistemo il primo corpo.
  if (root==NULL){ 

    root=malloc(sizeof(nodo_t)); 
    if (root == NULL){ perror("malloc"); exit(EXIT_FAILURE); }
    
root->massa=m;
root->x=x;
root->y=y;
}

      //Dal secondo corpo in poi
else{ root= inserisci_ricors(m, x, y,root, 0, 0, _s);}

return root;
}
      









void free_albero(nodo_t * root){

//se il puntatore alla radice è NULL, esco
  if (root==NULL) return;
//altrimenti richiamo la funzione e libero la memoria a partire dalle sue foglie
  else{

    if(root->NE != NULL) free_albero(root->NE);
    
    if(root->SE != NULL) free_albero(root->SE);
    
    if(root->NW != NULL) free_albero(root->NW);
    
    if(root->SW != NULL) free_albero(root->SW);
    
    else free(root);

return;
  }
}


//parametri x, y, root equivalenti a quelli della funzione massadi, parametri x_c, y_c, k equivalenti a quelli della funzione "inserisci_ricors"

double massadi_ricors(double x, double y, nodo_t * root, double x_c, double y_c, double k){

  k=k/2;

//se mi trovo in una foglia e le coordinate coincidono con quelle ricercate, ritorno la massa di quella foglia, altrimenti ritorno 0.
  if(root->NE==NULL && root->SE ==NULL && root->NW==NULL && root->SW == NULL){
    if(root->x==x && root->y==y) return root->massa;
    else return 0;
  }
//invece se non mi trovo in una foglia, in base alle coordinate ricercate, mi sposto ricorsivamente nei vari sottoquadranti andando di volta in volta a considerare il sottoquadrante che contiene le coordinate a cui sono interessato
  else{
    if(x>x_c){
//in ogni chiamata ricorsiva della funzione massadi_ricors modifico le coordinate del centro aggiungengo o sottraendo opportunamente la variabile k che ad ogni chiamata viene dimezzata
    if(y>=y_c) return massadi_ricors(x, y, root->NE, x_c + k , y_c + k, k);
    else return massadi_ricors(x, y, root->SE, x_c +k, y_c -k, k);
  }

  else {	//if x < x_c
    if(y>=y_c) return massadi_ricors(x, y, root->NW, x_c - k , y_c + k, k);
    else return massadi_ricors(x, y, root->SW, x_c - k, y_c -k, k);
  }
  }
}





double massadi(double x, double y, nodo_t * root){

//se il puntatore alla radice principale è NULL, esco dalla funzione
if(root==NULL) return 0;


//altrimenti richiamo la mia funzione assegnando come coordinate del centro (0,0) e come raggio _s.
else return massadi_ricors(x, y, root, 0, 0, _s);

}



#define G 6.67*pow(10, -11)

void stima_forza_ricors(double m, double x, double y, double *fx, double *fy, double theta, nodo_t * root, double s){

  s=s/2;
  double m2=0, x2=0, y2=0, r=0, F=0;

  //se la massa "m" che sto considerando , sulla quale devo stimare la forza, non c'è nel sottoalbero che sto considerando
  if(massadi(x, y, root) == 0){

	//prendo come massa quella del nodo in questione che chiamo m2
    m2=root->massa;
    if(m2==0) return;

    x2=root->x;
    y2=root->y;
	// e ne calcolo la distanza dalla massa m
    r = sqrt( pow(x-x2, 2) + pow(y-y2, 2) );

    //se il corpo m2 è una foglia calcolo la forza di m2 su m

    if(root->NE==NULL && root->SE ==NULL && root->NW==NULL && root->SW == NULL){
      F=(G*m*m2)/pow(r, 2);
      *fx+=F*(x2-x)/r;
      *fy+=F*(y2-y)/r;
    }


    // invece se il corpo è un centro di massa ma viene rispettata la condizione di soglia, calcolo la forza.
    else if(s/r < theta){

      F=(G*m*m2)/pow(r, 2);
      *fx+=F*(x2-x)/r;
      *fy+=F*(y2-y)/r;
    }
    //se invece non viene rispettata la condizione di soglia, richiamo la funzione nei suoi figli
    else{	
      stima_forza_ricors(m, x, y, fx, fy, theta, root->NE, s);
      stima_forza_ricors(m, x, y, fx, fy, theta, root->SE, s);
      stima_forza_ricors(m, x, y, fx, fy, theta, root->NW, s);
      stima_forza_ricors(m, x, y, fx, fy, theta, root->SW, s);

    }


  }//chiude l'if principale (if massadi==0)


//se invece massadi non torna 0 e quindi c'è la massa m nel sottoalbero considerato
 else{


	//se sono in una foglia e mi trova la massa vuol dire che sono nel mio stesso corpo m quindi torno indietro
    if(root->NE==NULL && root->SE ==NULL && root->NW==NULL && root->SW == NULL) return;


	//altrimenti, dato che mi trovo in un centro di massa che però comprende anche m,

    else{  //modifico la massa m2 da considerare la sua posizione, escludendo da esso m2
      m2=(root->NE->massa + root->SE->massa + root->NW->massa + root->SW->massa - m);
      x2 = ((root->NE->massa)*(root->NE->x) + (root->SE->massa)*(root->SE->x) + (root->NW->massa)*(root->NW->x) + (root->SW->massa)*(root->SW->x) - m*x) / m2;
      y2 = ((root->NE->massa)*(root->NE->y) + (root->SE->massa)*(root->SE->y) + (root->NW->massa)*(root->NW->y) + (root->SW->massa)*(root->SW->y) -m*y) / m2;

	//calcolo la distanza di questa m2 da m
      r = sqrt( pow(x-x2, 2) + pow(y-y2, 2) );


	//se la distanza del centro di massa da m soddisfa la condizione di soglia, ne calcolo la forza
      if(s/r < theta){

	F=(G*m*m2)/pow(r, 2);
	*fx+=F*(x2-y)/r;
	*fy+=F*(y2-y)/r;
      }
	//altrimenti procedo nei suoi figli
      else{	
	stima_forza_ricors(m, x, y, fx, fy, theta, root->NE, s);
	stima_forza_ricors(m, x, y, fx, fy, theta, root->SE, s);
	stima_forza_ricors(m, x, y, fx, fy, theta, root->NW, s);
	stima_forza_ricors(m, x, y, fx, fy, theta, root->SW, s);
      }

    }//chiude else (prima di modificare le coordinate di m2)
 }//chiude else principale.
  return;
}





void stima_forza(double x, double y, double*fx, double *fy, double theta, nodo_t * root){
*fx=0, *fy=0;
  double m=massadi(x, y, root);
  if(m==0) return;

  else stima_forza_ricors(m, x, y, fx, fy, theta, root, 4*_s);
//come ultimo argomento della funzione stima_forza_ricors metto 4*_s in quanto quando chiamo quella funzione, viene dimezzata quella quantità e l'"s" usato durante la funzione è 2*_s che è la dimensione di s per la radice principale 
return;
}
