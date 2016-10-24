#include <stdio.h>
#include <string.h>
#include <errno.h>
#include <getopt.h>
#include <unistd.h>
#include "extraction_motif.h"
#include <stdbool.h>

//**************************************************************************************
// Aide du programme
void usage ()
{ 
	fputs(
	"\n"
	"**************************************************************** \n"
	"*Extraction de motifs communs à plusieurs séquences biologiques*\n"
	"****************************************************************\n"
	"\n"
	"Usage : ./extraction_motif -f file.fa -d d_Levenshtein -q quorum  \n"
	"\n"
	"Options :\n",stderr);
	PRINT_OPTION("f","filename_seq","FILE","Nom du fichier FASTA contenant les séquences à étudier");
	PRINT_OPTION("m","longueur_motif","INT","Longueur du motif commun à extraire");
	PRINT_OPTION("d","d_Levenshtein","INT","Le nombre maximal d'erreurs autorisées entre occurrence et motif (distance de Levenshtein)");
	PRINT_OPTION("q","quorum","FLOAT"," Le pourcentage minimum de séquences qui doivent présenter une occurrence du motif ");
	PRINT_OPTION("s","save","FILE"," Nom du fichier dans lequel l'utilisateur veut sauvegarder les résultats de l'extraction");
	PRINT_OPTION("h","help","FILE","Aide du programme ");
	fputs("\n",stderr);
}

//**************************************************************************************
// libération de la mémoire allouée dynamiquement pour toutes les occurrences

void free_occurrence(TPOccurrence tete_occurrence)
{
	TPOccurrence occurrence_suivante = tete_occurrence->occ_next;
	free(tete_occurrence);
	if (occurrence_suivante != NULL){free_occurrence(occurrence_suivante);}
}
//**************************************************************************************
// libération de la mémoire allouée dynamiquement pour toutes les séquences 

void free_Sequence(TPSequence tete_sequence)
{
	TPSequence sequence_suivante = tete_sequence->seq_next;
	free_occurrence(tete_sequence->tete_occurrence);
	free(tete_sequence);
	if (sequence_suivante != NULL){free_Sequence(sequence_suivante);}
}

//**************************************************************************************
// libération de la mémoire allouée dynamiquement pour tous les motifs

void free_Motif(TPMotif tete_motif)
{
	TPMotif motif_suivant = tete_motif->motif_next;
	free(tete_motif->Motif);
	free_Sequence(tete_motif->tete_sequence);
	free(tete_motif);
	if (motif_suivant != NULL){free_Motif(motif_suivant);}
}

//**************************************************************************************
// libération de la mémoire allouée dynamiquement des motifs ne respectant pas le quorum

void free_Motif_quorum(TPMotif motif)
{
	free(motif->Motif);
	free_Sequence(motif->tete_sequence);
	free(motif);
}

//**************************************************************************************
// procédure de sécurité, vérifiant que la mémoire a bien été allouée pour un pointeur donné

void* _safeMalloc(const char* filename,int line, size_t size)
	{
	void* ptr= malloc(size);
	if(ptr==NULL)
		{
		fprintf(stderr,"[%s:%d] Probleme allocation mémoire",filename,line);
		exit(EXIT_FAILURE);\
		}	
	return ptr;
	}

//**************************************************************************************
// procédure permettant l'initialisation d'une nouvelle occurrence

void initialisation_occurrence(TPOccurrence* occurrence_current, int pos, int num_seq, int ins, int subst, int del, char last )
{
	(*occurrence_current)->pos=pos;
	(*occurrence_current)->seq=num_seq;
	(*occurrence_current)->ins=ins;
	(*occurrence_current)->subst=subst;
	(*occurrence_current)->del=del;
	(*occurrence_current)->last=last;
	(*occurrence_current)->occ_next = NULL;
}


//**************************************************************************************
// fonction permettant de prédire si on peut étendre le motif d'un caractère en plus, si ce n'est pas le cas la fonction renvoie une valeur de -1 pour la position suivante.

int calcul_position_caractere_suivant(TPSequence p_sequence, TPOccurrence p_occurrence, char** adr_Liste_seq, int i)
{
	int j=0 ;
	int longueur_sequence = strlen(adr_Liste_seq[ p_sequence->num_seq - 1]);
	j = p_occurrence->pos + i + p_occurrence->ins - p_occurrence->del - 1; 
	if (j == longueur_sequence) j= -1;  // si on arrive à la fin de la séquence, j prend la valeur de -1. 

	return j ; 	
}

//**************************************************************************************
// procédure permettant l'ajout d'un motif en tête d'une liste de motifs
void ajout_tete_motif (TPMotif* adr_tete_motif, char* word)
{
	TPMotif motif_current = safeMalloc(sizeof(TMotif));
	motif_current->motif_next = NULL;
	motif_current->tete_sequence = NULL;
	motif_current->Motif = safeMalloc(sizeof(char));
	strcpy(motif_current->Motif,word);
	motif_current->nb_seq_avec_occurrence = 0;
	motif_current->motif_next = *adr_tete_motif;
	*adr_tete_motif = motif_current;
}
//**************************************************************************************
// procédure permettant l'ajout d'un séquence en tête d'une liste de séquences
void ajout_tete_sequence (TPSequence* adr_tete_sequence, int num_sequence)
{
	TPSequence sequence_current = safeMalloc(sizeof(TSequence));
	sequence_current->seq_next = NULL;
	sequence_current->tete_occurrence = NULL;
	sequence_current->num_seq = num_sequence;
	sequence_current->seq_next = *adr_tete_sequence;
	*adr_tete_sequence = sequence_current;

}




//**************************************************************************************
// procédure permettant l'ajout d'un occurrence en tête d'une liste d'occurrences
void ajout_tete_occurrence (TPOccurrence* adr_tete_occurrence, TPOccurrence new_occurrence)
{

	TPOccurrence occurrence_current = safeMalloc(sizeof(TOccurrence));
	occurrence_current->occ_next = NULL;

	occurrence_current->occ_next = *adr_tete_occurrence;
	*adr_tete_occurrence = occurrence_current;
	

}



//**************************************************************************************
//procédure permettant d'affecter une occurrence à une position donnée 
void affectation (TPMotif* tete_motif, TPMotif prec_motif, TPSequence prec_sequence, TPOccurrence prec_occurrence, TPOccurrence new_occurrence)
{
	if (prec_occurrence != NULL)
	{
		prec_occurrence->occ_next = new_occurrence; 
		return;
	}

	//prec_occurrence == NULL
	if (prec_sequence != NULL)
	{
		prec_sequence->seq_next->tete_occurrence = new_occurrence; 
		return;
	}

	//prec_sequence == NULL
	if(prec_motif != NULL)
	{
		prec_motif->motif_next->tete_sequence->tete_occurrence = new_occurrence; 
		return;
	}
	(*tete_motif)->tete_sequence->tete_occurrence = new_occurrence;
}



//**************************************************************************************
//fonction permettant de savoir si un motif existe
bool existe_motif(TPMotif tete_motif, char* word, TPMotif* adr_prec_motif)
{

	TPMotif p = safeMalloc(sizeof(TMotif));
	p= tete_motif;
	*adr_prec_motif = NULL;

	while((p != NULL) && strcmp(p->Motif,word) < 0){
		
		*adr_prec_motif = p;
		 
		p = p->motif_next;

	}

	return ((p != NULL) && (strcmp(p->Motif,word) == 0)); //vrai
}



//**************************************************************************************
//fonction permettant de savoir si une séquence existe
bool existe_sequence(TPSequence tete_seq, int num_sequence, TPSequence* adr_prec_sequence)
{
	TPSequence p = safeMalloc(sizeof(TSequence));
	p = tete_seq;
	*adr_prec_sequence = NULL;

	while((p != NULL) && (p->num_seq < num_sequence)){
		
		*adr_prec_sequence = p;
		 
		p = p->seq_next;

	}
	
	return ((p != NULL) && (p->num_seq == num_sequence)); //vrai
}


//**************************************************************************************
//procédure permettant de renvoyer le dernier élément de la liste dans un prec_occurrence
void parcours_occurrence(TPOccurrence tete_occ, TPOccurrence* adr_prec_occurrence)
{

	TPOccurrence p = safeMalloc(sizeof(TOccurrence));
	p = tete_occ;

	*adr_prec_occurrence = NULL;

	while (p != NULL){
		*adr_prec_occurrence = p;
		p = p->occ_next;

	}

}



//**************************************************************************************
// fonction permettant de vérifier l'existence d'un motif, d'une séquence et de renvoyer la position à la fin de la liste d'occurrence si la séquence existe
bool existe(TPMotif tete_motif, TPOccurrence new_occurrence, char* word, int num_sequence, int pos_occurrence, bool *adr_presence_motif, bool *adr_presence_sequence, bool *adr_presence_occurrence, TPMotif* adr_prec_motif, TPSequence* adr_prec_sequence, TPOccurrence* adr_prec_occurrence)
{
	TPMotif motif_current = safeMalloc(sizeof(TMotif));
	motif_current = NULL;
	*adr_presence_motif = false;
	*adr_presence_sequence = false;
	
	if(!existe_motif(tete_motif, word, adr_prec_motif))
	{

		return false;
	}

	*adr_presence_motif = true;

	if (*adr_prec_motif == NULL)
	{

		motif_current = tete_motif;
	}
	else{
		motif_current = (*adr_prec_motif)->motif_next;

	}

	if((motif_current != NULL) && existe_sequence(motif_current->tete_sequence, num_sequence, adr_prec_sequence))
	{

		*adr_presence_sequence = true;
		if(*adr_prec_sequence == NULL)
		{
			parcours_occurrence(motif_current->tete_sequence->tete_occurrence, adr_prec_occurrence);
	
		}
			 

		else if(*adr_prec_sequence != NULL) 
		{	
			
			parcours_occurrence((*adr_prec_sequence)->seq_next->tete_occurrence, adr_prec_occurrence);
	
		}


	}
	return false;
}




//**************************************************************************************
//procédure permettant d'ajouter un motif dans la liste de motifs
void creation_motif (TPMotif* adr_tete_motif, char* word, TPMotif prec_motif)
{

	if (prec_motif == NULL)
	{
		ajout_tete_motif(adr_tete_motif, word);
	}
	else{
		ajout_tete_motif(&prec_motif->motif_next, word);
	}
}



//**************************************************************************************
//procédure permettant d'ajouter une séquence dans une liste de séquences
void creation_ligne_sequence (TPSequence* adr_tete_sequence, int num_sequence, TPSequence prec_sequence)
{
	if (prec_sequence == NULL)
	{
		ajout_tete_sequence(adr_tete_sequence, num_sequence);
	}
	else{
		ajout_tete_sequence(&prec_sequence->seq_next, num_sequence);
	}
}





//**************************************************************************************
//procédure permettant d'ajouter une occurrence dans une liste d'occurrences
void creation_col_occurrence(TPOccurrence* adr_tete_occurrence, TPOccurrence new_occurrence, TPOccurrence prec_occurrence)
{

	if (prec_occurrence == NULL)
	{
		ajout_tete_occurrence(adr_tete_occurrence, new_occurrence);

	}
	else{
		ajout_tete_occurrence(&prec_occurrence->occ_next, new_occurrence);
	}
}





//**************************************************************************************
//procédure permettant de de vérifier la présence de motif et de séquence, si elle existe récupérer l'adresse respectivement du motif précédent et de la séquence précédente, et d'ajouter la nouvelle occurrence
void creer_col_occurrence(TPMotif* tete_motif, char* word, int num_sequence, bool presence_motif, bool presence_sequence, TPMotif prec_motif, TPSequence prec_sequence, TPOccurrence prec_occurrence, TPOccurrence new_occurrence)
{

	TPMotif motif_current = safeMalloc(sizeof(TMotif));

	if (!presence_motif)
	{

		creation_motif(tete_motif, word, prec_motif);
		prec_sequence = NULL;
	}

	if (prec_motif == NULL)
	{	
		motif_current = *tete_motif;

	}
	else{

		motif_current = prec_motif->motif_next;
	}

	if (!presence_sequence)
	{
		creation_ligne_sequence(&motif_current->tete_sequence, num_sequence, prec_sequence);
		prec_occurrence = NULL;
		motif_current->nb_seq_avec_occurrence++;

	}

	if(prec_sequence == NULL){

		creation_col_occurrence(&motif_current->tete_sequence->tete_occurrence, new_occurrence, prec_occurrence);


	}else{
		creation_col_occurrence(&prec_sequence->seq_next->tete_occurrence, new_occurrence, prec_occurrence);
	}

	affectation(tete_motif, prec_motif, prec_sequence, prec_occurrence, new_occurrence);


}




//**************************************************************************************
//procédure permettant de placer une nouvelle occurrence dans le dictionnaire de motif
void set(TPMotif* tete_motif, char* word, int num_sequence, int pos_occurrence, TPOccurrence new_occurrence)
{
	TPMotif prec_motif= safeMalloc(sizeof(TMotif));
	prec_motif = NULL;
	TPSequence prec_sequence=safeMalloc(sizeof(TSequence));
	prec_sequence = NULL;
	TPOccurrence prec_occurrence = safeMalloc(sizeof(TOccurrence));
	prec_occurrence = NULL;
	bool presence_motif;
	bool presence_sequence;
	bool presence_occurrence;

	bool valeur_remplacee_est_val_def = !existe(*tete_motif, new_occurrence, word, num_sequence, pos_occurrence, &presence_motif, &presence_sequence, &presence_occurrence, &prec_motif, &prec_sequence, &prec_occurrence);

	if(valeur_remplacee_est_val_def){

		creer_col_occurrence(tete_motif, word, num_sequence, presence_motif, presence_sequence, prec_motif, prec_sequence, prec_occurrence, new_occurrence);
	}
	else{
		affectation(tete_motif, prec_motif, prec_sequence, prec_occurrence, new_occurrence);
	}

}




//**************************************************************************************
//fonction permettant de compter le nombre de séquences du fichier donné en entrée. 

int compte_nb_sequences(char* filename_seq) 
{
	FILE* file_seq=NULL; 
	int nombre_sequences=0;
	char seq[FILENAME_MAX];  

	if (filename_seq != NULL) {
		file_seq=fopen(filename_seq,"r");
	}	

	if ( file_seq == NULL) { 
		fprintf(stderr, "Impossible d'ouvrir le fichier \n"); EXIT_FAILURE;
	}
	
	while(fgets(seq,FILENAME_MAX, file_seq ) != NULL) {
		if((seq[0] == '>') || (seq[0] == '\n'))continue;  // ne compte pas les noms des séquences et les lignes vides. 
		nombre_sequences++;	
	}


	fclose(file_seq);

	return nombre_sequences;
}




//**************************************************************************************
// procédure permettant d'ouvrir le fichier fasta donné en argument et de stocker les séquences qu'il contient dans un tableau dynamique ( Liste_seq ).

void OpenFile(char* filename_seq, char** Liste_seq) 
{
	FILE* file_seq=NULL; 
	char seq[FILENAME_MAX];
	int cpt_seq =0 ;

	if (filename_seq != NULL) {
		file_seq=fopen(filename_seq,"r");
	}	

	if ( file_seq == NULL) {
			fprintf(stderr, "Cannot open file \n");	
	}

	while(fgets(seq,FILENAME_MAX, file_seq ) != NULL) // lecture du fichier ligne par ligne. 
   	{ 	 
		size_t line_len=strlen(seq);
		if(line_len==0) continue;	// ne prend pas en compte les lignes vides
		if((seq[0] == '>') || (seq[0] == '\n'))continue;  // ne prend pas en compte les sauts à la ligne et les lignes avec le nom des séquences. 
		if(seq[line_len - 1] == '\n') { seq[line_len - 1] = '\0';}   // transforme le saut à la ligne en \0.	
		else { fprintf(stderr,"Cannot read file [%s]\n",strerror(errno));}
		
		Liste_seq[cpt_seq] = (char*)malloc(sizeof(char)*(line_len+1)); // allocation de la mémoire pour la nouvelle séquence
		strcpy(Liste_seq[cpt_seq], seq);  // stockage de la nouvelle séquence dans Liste_seq. 

		cpt_seq++;	
    	}

	fclose(file_seq);	

}




//**************************************************************************************
// procédure permettant de créer le dictionnaire de motif vide
void creer_dictionnaire__motif_vide(TPMotif* tete_motif, char* word, int nombre_sequences, char** Liste_seq, int d_Levenshtein)
{
	int i =0, nb_max_insertions = 0;
	for (i=0; i < nombre_sequences; i++)
	{
		int len = strlen(Liste_seq[i]);
		int j = 0;
		for (j=0; j<len;j++)
		{
			TPOccurrence occurrence_current = safeMalloc(sizeof(TOccurrence));
			initialisation_occurrence(&occurrence_current, j+1, i+1, 0, 0, 0, 'm' );
			set(tete_motif, word, i+1, j+1, occurrence_current);
			
			nb_max_insertions = min(d_Levenshtein, (len-j));
			int r=1;
			for (r=1; r<=nb_max_insertions; r++)
			{
				TPOccurrence occurrence_current = safeMalloc(sizeof(TOccurrence));
				initialisation_occurrence(&occurrence_current, j+1, i+1, r, 0, 0, 'i' );
				set(tete_motif, word, i+1, j+1, occurrence_current);
			}

		}
	}

}



//**************************************************************************************
// procédure permettant de supprimer la tête d'une liste de motifs
void supprimer_tete_motif(TPMotif* adr_tete)
	{
		TPMotif ancienne_tete = *adr_tete;
		*adr_tete = ancienne_tete -> motif_next;
		free_Motif_quorum(ancienne_tete);

	}




//**************************************************************************************
// fonction permettant l'élongation du motif
char* elongation_motif(TMotif* p_motif, char nucleotide_ajout, int i)
{
	char* motif_alonge = NULL;
	if (strcmp(p_motif->Motif, "") == 0)
		{
			motif_alonge = safeMalloc(sizeof(char)*2);
		}
		else
		{
			motif_alonge = safeMalloc(sizeof(char)*(i+1));
			strcpy(motif_alonge,p_motif->Motif);
		}
		char* new_nucleotide = safeMalloc(sizeof(char));
		new_nucleotide[0] = nucleotide_ajout;
		new_nucleotide[1] = '\0';
		strcat(motif_alonge, new_nucleotide);
		free(new_nucleotide);

		return motif_alonge;
}





// ENREGISTREMENT DU RESULTATS DE L'EXTRACTION DE MOTIFS DANS UNE FICHIER
// si un nom de fichier est donnée en argument, on crée ce fichier
void exportation_fichier(char* nom_fichier_sauvegarde, TPMotif new_tete_motif, char* choixMotif, int quorum, int nombre_sequences)
{
	TPMotif parcours_motif = safeMalloc(sizeof(TMotif));
	TPSequence t_sequence = safeMalloc(sizeof(TMotif));
	TPOccurrence t_occurrence = safeMalloc(sizeof(TMotif));

	FILE* fichier = NULL;
	fichier = fopen(nom_fichier_sauvegarde, "a"); // écrire dans le fichier en mode "append"

	if (fichier != NULL){
		parcours_motif = new_tete_motif;
		if (new_tete_motif == NULL){fprintf(stderr, "Le dictionnaire est vide.\n");}
		while(parcours_motif != NULL)
		{
			fprintf(fichier, "--------------------------- %s --------------------------------------\n", parcours_motif->Motif);
			t_sequence = parcours_motif->tete_sequence;
			while(t_sequence != NULL)
			{
				fprintf(fichier,"quorum : %f\n", (float)parcours_motif->nb_seq_avec_occurrence /  (float)nombre_sequences);
				fprintf(fichier, "++++ seq %d ++++\n", t_sequence->num_seq);
				t_occurrence = t_sequence->tete_occurrence;
				while(t_occurrence != NULL)
				{
					fprintf(fichier, "pos %d, ins %d, del %d, subst %d, last %c\n", t_occurrence->pos,  t_occurrence->ins,  t_occurrence->del, t_occurrence->subst, t_occurrence->last);
					t_occurrence = t_occurrence->occ_next;
				}
				t_sequence = t_sequence->seq_next;
			}
			parcours_motif=parcours_motif->motif_next;
		}
	}
	else  fprintf(fichier, "impossible de sauvegarder le resultat de l'extraction dans un fichier");
	fclose(fichier);
}





//**************************************************************************************
// fonction permettant de choisir un motif à consulter
void choisir_motif(char* choixMotif, int long_motif, TPMotif* new_tete_motif, int nombre_sequences)
{
	TPMotif parcours_motif = safeMalloc(sizeof(TMotif));
	TPSequence t_sequence = safeMalloc(sizeof(TMotif));
	TPOccurrence t_occurrence = safeMalloc(sizeof(TMotif));

	fprintf(stderr, "choisissez un motif à consulter:\n");
	scanf("%s", choixMotif);
	fprintf(stderr, "taille motif: %d\n", long_motif);
	
	parcours_motif = *new_tete_motif;
	if (*new_tete_motif == NULL){fprintf(stderr, "Le dictionnaire est vide.\n");return;}
		while((parcours_motif != NULL) && (strcmp(parcours_motif->Motif,choixMotif)!=0))
		{
			parcours_motif=parcours_motif->motif_next;
		}

	if (parcours_motif == NULL){fprintf(stderr, "Le motif est absent du dictionnaire.\n");return;}
	fprintf(stderr, "--------------------------- %s --------------------------------------\n", parcours_motif->Motif);
	t_sequence = parcours_motif->tete_sequence;
	while(t_sequence != NULL)
	{

		fprintf(stderr,"quorum : %f\n", (float)parcours_motif->nb_seq_avec_occurrence /  (float)nombre_sequences);
		
		fprintf(stderr, "++++ seq %d ++++\n", t_sequence->num_seq);
		
		t_occurrence = t_sequence->tete_occurrence;
		while(t_occurrence != NULL)
		{
			fprintf(stderr, "pos %d, ins %d, del %d, subst %d, last %c\n", t_occurrence->pos,  t_occurrence->ins,  t_occurrence->del, t_occurrence->subst, t_occurrence->last);
			t_occurrence = t_occurrence->occ_next;
		}
		t_sequence = t_sequence->seq_next;
	}
}





//**************************************************************************************
// Algorithme générale pour l'extraction de motifs communs à plusieurs séquences

void recherche_de_motifs(TPMotif* new_tete_motif, int nombre_sequences, char* filename_seq, int long_motif, float quorum, int d_Levenshtein)
{

	//*********************************************************
	// Récupération des séquences dans le fichier FASTA
	
	char* Liste_seq[nombre_sequences-1];
	char** p_Liste_seq = &Liste_seq[0];
	int i;

	OpenFile(filename_seq, p_Liste_seq);




	//********************************************************
	// Dictionnaire vide
	

	TPMotif motif_vide = safeMalloc(sizeof(TMotif));
	motif_vide = NULL;
	
	creer_dictionnaire__motif_vide(&motif_vide, "", nombre_sequences, Liste_seq, d_Levenshtein);
	

	//*********************************************************
	// Création d'un tableau contenant les 4 différents acides nucléiques avec lesquels on peut allonger le motif. 

	char ensemble_nucleotides[3];
	ensemble_nucleotides[0] = 'A';
	ensemble_nucleotides[1] = 'C';
	ensemble_nucleotides[2] = 'G';
	ensemble_nucleotides[3] = 'T';


	//*********************************************************
	// Début de l'algorithme : Extraction de motifs communs à plusieurs séquences biologiques

	*new_tete_motif = motif_vide;

	for (i = 0; i <= long_motif-1; ++i)
	{
		TPMotif p_motif = safeMalloc(sizeof(TMotif));
		TPMotif tete_new_dict = safeMalloc(sizeof(TMotif));
		tete_new_dict = NULL;
		p_motif = *new_tete_motif;

		//creation dictionnaire de taille i:

		TPSequence p_sequence = safeMalloc(sizeof(TSequence));
		TPOccurrence p_occurrence = safeMalloc(sizeof(TOccurrence));
		
		char* newmotif = NULL;
		char* new_nucleotide = NULL;


		while (p_motif != NULL) //pour chaque motif M
		{	
			p_sequence = p_motif->tete_sequence;
			
			while (p_sequence != NULL) //pour chaque sequence
			{
				p_occurrence = p_sequence->tete_occurrence;

				while (p_occurrence != NULL) //pour chaque occurrence
				{	
					//calcul de la position suivante
					int j = calcul_position_caractere_suivant(p_sequence, p_occurrence, p_Liste_seq, i); //calcul de la position j du caractère suivant	
					

					//cas de la substitution?
					switch(p_occurrence->last)
					{
						case 'm':
							if ((j!= -1) && (p_occurrence->subst+p_occurrence->ins+p_occurrence->del < d_Levenshtein))
							{
								int k; 
								for (k = 0; k < 4; ++k)
								{
									if (ensemble_nucleotides[k] != p_Liste_seq[p_sequence->num_seq-1][j])
									{
										//élongation du motif avec les ensembles de nucléotides
										newmotif = elongation_motif(p_motif, ensemble_nucleotides[k], i);

										TPOccurrence occurrence_current = safeMalloc(sizeof(TOccurrence));
										//initialisation de la nouvelle occurrence
										initialisation_occurrence(&occurrence_current, p_occurrence->pos , p_sequence->num_seq, p_occurrence->ins, p_occurrence->subst+1,  p_occurrence->del, 's');
										//enregistrement de la nouvelle occurrence dans le nouveau dictionnaire
										set(&tete_new_dict, newmotif, p_sequence->num_seq, p_occurrence->pos, occurrence_current);
									}
								}
							}
							break;
							
						case 'i':
						case 's':
							if ((j!= -1) && (p_occurrence->subst+p_occurrence->ins+p_occurrence->del < d_Levenshtein))
							{
								int k;
								for (k = 0; k < 4; ++k)
								{
									if ((ensemble_nucleotides[k] != p_Liste_seq[p_sequence->num_seq-1][j]) && (ensemble_nucleotides[i] != p_Liste_seq[p_sequence->num_seq-1][j-1]))
									{
										//élongation du motif avec les ensembles de nucléotides
										newmotif = elongation_motif(p_motif, ensemble_nucleotides[k], i);

										TPOccurrence occurrence_current = safeMalloc(sizeof(TOccurrence));
										//initialisation de la nouvelle occurrence
										initialisation_occurrence(&occurrence_current, p_occurrence->pos , p_sequence->num_seq, p_occurrence->ins, p_occurrence->subst+1,  p_occurrence->del, 's');
										//enregistrement de la nouvelle occurrence dans le nouveau dictionnaire
										set(&tete_new_dict, newmotif, p_sequence->num_seq, p_occurrence->pos, occurrence_current);
									}
								}
							}
							break;
					}






					//cas de la deletion?
					switch(p_occurrence->last)
					{
						case 'm':
						case 's':
							if (p_occurrence->subst+p_occurrence->ins+p_occurrence->del < d_Levenshtein)
							{
								int k; 
								for (k = 0; k < 4; ++k)
								{
									if (ensemble_nucleotides[k] != p_Liste_seq[p_sequence->num_seq-1][j-1])
									{
										//élongation du motif avec les ensembles de nucléotides
										newmotif = elongation_motif(p_motif, ensemble_nucleotides[k], i);

										TPOccurrence occurrence_current = safeMalloc(sizeof(TOccurrence));
										//initialisation de la nouvelle occurrence
										initialisation_occurrence(&occurrence_current, p_occurrence->pos , p_sequence->num_seq, p_occurrence->ins, p_occurrence->subst,  p_occurrence->del+1, 'd');
										//enregistrement de la nouvelle occurrence dans le nouveau dictionnaire
										set(&tete_new_dict, newmotif, p_sequence->num_seq, p_occurrence->pos, occurrence_current);
									
									}
								}
							}

						break;

						case 'd':
							if (p_occurrence->subst+p_occurrence->ins+p_occurrence->del < d_Levenshtein)
							{
								int k;
								for (k = 0; k < 4; ++k)
								{
									if (ensemble_nucleotides[k] != p_Liste_seq[p_sequence->num_seq-1][j])
									{
										//élongation du motif avec les ensembles de nucléotides
										newmotif = elongation_motif(p_motif, ensemble_nucleotides[k], i);

										TPOccurrence occurrence_current = safeMalloc(sizeof(TOccurrence));
										//initialisation de la nouvelle occurrence
										initialisation_occurrence(&occurrence_current, p_occurrence->pos , p_sequence->num_seq, p_occurrence->ins, p_occurrence->subst,  p_occurrence->del+1, 'd');
										//enregistrement de la nouvelle occurrence dans le nouveau dictionnaire
										set(&tete_new_dict, newmotif, p_sequence->num_seq, p_occurrence->pos, occurrence_current);
									
									}
								}
							}

						break;
					}
					


					// cas du match ?
					bool match_possible = true;
					if (j == -1) {
						match_possible = false;
						p_occurrence = p_occurrence->occ_next;
						continue; 
					}

					if ((p_occurrence->last == 'i') && (p_Liste_seq[p_sequence->num_seq-1][j] == p_Liste_seq[p_sequence->num_seq-1][j-1])){p_occurrence = p_occurrence->occ_next;continue;}

					//élongation du motif
					newmotif = elongation_motif(p_motif, p_Liste_seq[p_sequence->num_seq-1][j], i);

					


					//fin élongation de motif
					TPOccurrence occurrence_current = safeMalloc(sizeof(TOccurrence));
					//initialisation de la nouvelle occurrence
					initialisation_occurrence(&occurrence_current, p_occurrence->pos , p_sequence->num_seq, p_occurrence->ins, p_occurrence->subst,  p_occurrence->del, 'm');
					//enregistrement de la nouvelle occurrence dans le nouveau dictionnaire
					set(&tete_new_dict, newmotif, p_sequence->num_seq, p_occurrence->pos, occurrence_current);



					

					//cas d'insertions ?

					int nb_max_insertions_theoriques = 0;
					int nb_caractere_dispo_a_droite = 0;
					int nb_max_insertions=0;
					int r=0;

					nb_max_insertions_theoriques = d_Levenshtein - (p_occurrence->subst+p_occurrence->ins+p_occurrence->del);	
					if (j == -1) continue;
					nb_caractere_dispo_a_droite = strlen(p_Liste_seq[p_sequence->num_seq-1]) - j - 1;
					nb_max_insertions = min(nb_caractere_dispo_a_droite,nb_max_insertions_theoriques);
					
					//élongation du motif
					newmotif = elongation_motif(p_motif, p_Liste_seq[p_sequence->num_seq-1][j], i);

					for (r=1; r<=nb_max_insertions; r++){
						
						TPOccurrence occurrence_current = safeMalloc(sizeof(TOccurrence));
						//initialisation de la nouvelle occurrence
						initialisation_occurrence(&occurrence_current, p_occurrence->pos , p_sequence->num_seq, p_occurrence->ins +1, p_occurrence->subst,  p_occurrence->del, 'i');
						//enregistrement de la nouvelle occurrence dans le nouveau dictionnaire	
						set(&tete_new_dict, newmotif, p_sequence->num_seq, p_occurrence->pos, occurrence_current);
						j++;
					}
	
					p_occurrence = p_occurrence->occ_next;
					free(new_nucleotide);	
				}
				p_sequence = p_sequence->seq_next;
			}
			p_motif = p_motif->motif_next;

		}
		



		// libération de la mémoire allouée dynamiquement pour le dictionnaire i-1

		if (*new_tete_motif != NULL)
		{
			free_Motif(*new_tete_motif);
		}
		
		*new_tete_motif = tete_new_dict;




		// Vérification si les motifs vérifient le quorum. 

		TPMotif motif_verification_quorum = safeMalloc(sizeof(TMotif));
		TPMotif	motifprec_verification_quorum = safeMalloc(sizeof(TMotif));
		motif_verification_quorum =  *new_tete_motif;
		motifprec_verification_quorum = NULL;
		while(motif_verification_quorum != NULL)
		{
			if (((float)motif_verification_quorum->nb_seq_avec_occurrence /  (float)nombre_sequences) < quorum){
				if (motifprec_verification_quorum != NULL){
					supprimer_tete_motif(&motifprec_verification_quorum->motif_next);
					motif_verification_quorum = motifprec_verification_quorum->motif_next;
				}
				else {
					supprimer_tete_motif(new_tete_motif);
					motif_verification_quorum = *new_tete_motif;
				}		

			}	
			else {
				motifprec_verification_quorum = motif_verification_quorum;
				motif_verification_quorum =motif_verification_quorum ->motif_next;	
			}	
		}
		free(motif_verification_quorum);
	}
}




//**************************************************************************************
//  Affichage du motif choisi par l'utilisateur ainsi que la liste des occurrences de ce motif

void affichage_liste_motifs(TPMotif* new_tete_motif, int long_motif)
{
	TPMotif parcours_motif = safeMalloc(sizeof(TMotif));
	parcours_motif = *new_tete_motif;
	if (new_tete_motif == NULL){ fprintf(stderr, "Le dictionnaire est vide.\n"); return;}
	fprintf(stderr, "------------- motifs de taille %d -------------\n", long_motif);
	
	if(parcours_motif == NULL){fprintf(stderr, "\nLe dictionnaire est vide\n");}
	while(parcours_motif != NULL)
	{
		fprintf(stderr, "%s\n", parcours_motif->Motif);
		parcours_motif=parcours_motif->motif_next;
	}
	free(parcours_motif);
}

