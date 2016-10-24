#include <stdio.h>
#include <string.h>
#include <errno.h>
#include <getopt.h>
#include <unistd.h>
#include "extraction_motif.h"
#include <stdbool.h>


// Début du programme
int main(int argc, char** argv)
{

	char* filename_seq = safeMalloc(sizeof(char));
	char* nom_fichier_sauvegarde = safeMalloc(sizeof(char));
	filename_seq = NULL;
	nom_fichier_sauvegarde = NULL;
	int long_motif, d_Levenshtein = 0;
	float quorum = 0;
	int nombre_sequences;


	//*********************************************************
	// getopt_long permettant de récuperer les options entrées en arguments du programme. 

	int c=0, option_index = 0;

        static struct option long_options[] = {
       	    {"help", no_argument, NULL, 'h'},
            {"filename_seq",  required_argument, 0, 'f'},
            {"longueur_motif",  required_argument, 0, 'm'},
            {"d_Levenshtein",  required_argument, 0, 'd'},
            {"quorum",  required_argument, 0, 'q'},
	    {"save",  required_argument, 0, 's'},
            {0,		0,	  0,	 0 }
        };

	while ((c = getopt_long(argc, argv,"hf:m:d:q:s:", long_options, &option_index)) != -1) {

		switch(c) {

			case 'h': usage();return EXIT_SUCCESS; break;
			
			case 'f': filename_seq = optarg; break;
			 
			case 'm': long_motif = atoi(optarg); break;
			
			case 'd': d_Levenshtein = atoi(optarg); break;

			case 'q': quorum =(float)atof(optarg); break;
			
			case 's': nom_fichier_sauvegarde = optarg; break;
			 
			case '?': fprintf(stderr, "ERROR: option -%c is undefined\n", optopt); return EXIT_FAILURE; break;
			
			case ':': fprintf(stderr, "ERROR: option -%c requires an argument\n",optopt); return EXIT_FAILURE; break;
			 	
			default : usage(); return EXIT_FAILURE; break;
	
		}
	}
	
	if (filename_seq == NULL) {  fprintf(stderr, "ERROR: Aucun fichier entrée en paramètre \n"); return EXIT_FAILURE;}

	TPMotif* new_tete_motif = safeMalloc(sizeof(TMotif));
	nombre_sequences = compte_nb_sequences(filename_seq);
	//lancement de la procédure de motifs
	recherche_de_motifs(new_tete_motif, nombre_sequences, filename_seq, long_motif, quorum, d_Levenshtein);



	//menu affiché à l'utilisateur
	int choixDuMenu;
	char choixMotif[long_motif+1];

	while(choixDuMenu != 3)
	{
		fprintf(stderr, "\n\n\t********************************************************\n");
		fprintf(stderr, "\t***********************   MENU   ***********************\n");
		fprintf(stderr, "\t********************************************************\n\n");
		fprintf(stderr, "\t1. afficher la liste des motifs\n");
		fprintf(stderr, "\t2. choisir un motif et consulter a liste des occurrences de ce motif\n");
		fprintf(stderr, "\t3. sortir du programme\n");
		
		scanf("%d", &choixDuMenu);

		switch(choixDuMenu){

			case 1:
				affichage_liste_motifs(new_tete_motif, long_motif);
				break;

			case 2:
				choisir_motif(choixMotif, long_motif, new_tete_motif, nombre_sequences);
				if (nom_fichier_sauvegarde != NULL) {
					exportation_fichier(nom_fichier_sauvegarde, *new_tete_motif, choixMotif, quorum, nombre_sequences);
				}
				break;
			default:
				return EXIT_FAILURE; break;
				
		}	
	}


	// liberation de la mémoire pour le dernier dictionnaire. 
	if (*new_tete_motif != NULL){free_Motif(*new_tete_motif);}

} // fin du programme.
