#include <stdio.h>
#include <stdlib.h>
#include <stdbool.h>

typedef struct TOccurrence TOccurrence;
struct TOccurrence  {
	int pos;
	int seq;
	int ins;
	int subst;
	int del;
	char last;
	struct TOccurrence *occ_next;
};


typedef struct TSequence TSequence;
struct TSequence{
	int num_seq;
	struct TSequence *seq_next;
	struct TOccurrence *tete_occurrence;
} ;


typedef struct TMotif TMotif;
struct TMotif{
	char* Motif;
	int nb_seq_avec_occurrence;
	struct TSequence *tete_sequence;
	struct TMotif *motif_next;
} ;

typedef TMotif* TPMotif;
typedef TSequence* TPSequence;
typedef TOccurrence* TPOccurrence;


#define safeMalloc(S) _safeMalloc(__FILE__,__LINE__,S)

//calcul du minimum entre deux entiers : 
#define min(a,b) (a<=b?a:b)

// Macro afin de mettre en forme les différentes options dans l'aide (-h). 
#define PRINT_OPTION(SHORT,LONG,TYPE,DEF) \
	fputs("   -" SHORT ",--" LONG "   " TYPE "  " DEF "\n",stderr)   




//*****************************************************************************************
// Prototypes de toutes les fonctions utilisées : 

void usage ();
void free_occurrence(TPOccurrence tete_occurrence);
void free_Sequence(TPSequence tete_sequence);
void free_Motif(TPMotif tete_motif);
void free_Motif_quorum(TPMotif motif);
void* _safeMalloc(const char* filename,int line, size_t size);
void initialisation_occurrence(TPOccurrence* occurrence_current, int pos, int num_seq, int ins, int subst, int del, char last );
int calcul_position_caractere_suivant(TPSequence p_sequence, TPOccurrence p_occurrence, char** adr_Liste_seq, int i);
void ajout_tete_motif (TPMotif* adr_tete_motif, char* word);
void ajout_tete_sequence (TPSequence* adr_tete_sequence, int num_sequence);
void ajout_tete_occurrence (TPOccurrence* adr_tete_occurrence, TPOccurrence new_occurrence);
void affectation (TPMotif* tete_motif, TPMotif prec_motif, TPSequence prec_sequence, TPOccurrence prec_occurrence, TPOccurrence new_occurrence);
bool existe_motif(TPMotif tete_motif, char* word, TPMotif* adr_prec_motif);
bool existe_sequence(TPSequence tete_seq, int num_sequence, TPSequence* adr_prec_sequence);
void existe_occurrence(TPOccurrence tete_occ, TPOccurrence* adr_prec_occurrence);
bool existe(TPMotif tete_motif, TPOccurrence new_occurrence, char* word, int num_sequence, int pos_occurrence, bool *adr_presence_motif, bool *adr_presence_sequence, bool *adr_presence_occurrence, TPMotif* adr_prec_motif, TPSequence* adr_prec_sequence, TPOccurrence* adr_prec_occurrence);
void creation_motif (TPMotif* adr_tete_motif, char* word, TPMotif prec_motif);
void creation_ligne_sequence (TPSequence* adr_tete_sequence, int num_sequence, TPSequence prec_sequence);
void creation_col_occurrence(TPOccurrence* adr_tete_occurrence, TPOccurrence new_occurrence, TPOccurrence prec_occurrence);
void creer_col_occurrence(TPMotif* tete_motif, char* word, int num_sequence, bool presence_motif, bool presence_sequence, TPMotif prec_motif, TPSequence prec_sequence, TPOccurrence prec_occurrence, TPOccurrence new_occurrence);
void set(TPMotif* tete_motif, char* word, int num_sequence, int pos_occurrence, TPOccurrence new_occurrence);
int compte_nb_sequences(char* filename_seq);
void OpenFile(char* filename_seq, char** Liste_seq); 
void creer_dictionnaire__motif_vide(TPMotif* tete_motif, char* word, int nombre_sequences, char** Liste_seq, int d_Levenshtein);
void supprimer_tete_motif(TPMotif* adr_tete);
char* elongation_motif(TMotif* p_motif, char nucleotide_ajout, int i);
void recherche_de_motifs(TPMotif* new_tete_motif, int nombre_sequences, char* filename_seq, int long_motif, float quorum, int d_Levenshtein);
void affichage_liste_motifs(TPMotif* new_tete_motif, int long_motif);
void exportation_fichier(char* nom_fichier_sauvegarde, TMotif* new_tete_motif, char* choixMotif, int quorum, int nombre_sequences);
void choisir_motif(char* choixMotif, int long_motif, TPMotif* new_tete_motif, int nombre_sequences);
