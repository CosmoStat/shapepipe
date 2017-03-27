
#ifndef _phil
#define _phil

#include <iostream>

#define _mac_assertmsg(TST,MSG) \
        ((TST) ? (void)0        \
	       : (cerr << __FILE__ << "(" << __LINE__ << \
                ") " << endl << "Assertion failed: " #TST << \
                 MSG << endl,abort()))

#define _mac_invariant(TST,MSG) \
        _mac_assertmsg (TST, " - Invariant : " << MSG)
#define _mac_precondition(TST,MSG) \
        _mac_assertmsg (TST, " - Precondition : " << MSG)
#define _mac_postcondition(TST,MSG) \
        _mac_assertmsg (TST, " - Postcondition : " << MSG)
	
	
	
	
/*******************************************************************************
Graph call
*******************************************************************************/	

#include <iostream>
#include <fstream>
#include <string>
#include <stdlib.h>

using std::string;
using std::ostream;
using std::cout;
using std::ofstream;
using std::endl;

// Pour plus de facilité d'utilisation :
#define TRACEUR	Traceur::GetSingleton()
#define INDENT	Indent::GetSingIdent()


class Indent {
private:
   int             _Decal;       // Current Indentation value
   static Indent*  _SingIdent;   // Unique Objet de la classe
public:
   int  getIndent() {return _Decal;}
   void incIndent() {_Decal+=2;}
   void incIndent(int Dec) {_Decal+=Dec;}
   void decIndent() {_Decal-=2;}
   void decIndent(int Dec) {_Decal-=Dec;}
   static Indent& GetSingIdent() {       // Accesseur statique
      if (_SingIdent != NULL) return *_SingIdent;
      _SingIdent = new Indent();
      _SingIdent->_Decal=0;
      return *_SingIdent;
   };	   
};

class Traceur;

class Trace {
	
	friend class Traceur;

public:
	// Type d'information à Tracer :
	enum TYPE {
		ERROR,		// Erreur grave.
		EXCEPTION,	// Exception, parfois, ca n'est pas une erreur.
		METIER,		// Information métier à tracer.
		MEMORY, 	// Construction et destruction d'objets.
		TECHNIQUE	// Info technique : Entrée et sortie de méthode...
	};
 
	// Niveau d'importance de la trace :
	enum NIVEAU {
		IMPORTANT,	// L'information doit oblig être relatée !
		VERBOSE,	// Information relatée en mode verbose
		DEBUG           // Info affichée en mode debug
	};
	
private:
// Attributs :
	TYPE	Type;		// Type d'information à Tracer.
	string	NomClasse;	// Nom de la classe de l'objet.
	string	NomMethode;	// Méthode ou la trace est crée.
	string	Message;	// Message de la trace.
	void*	Objet;		// Objet expéditeur de la trace.	
	NIVEAU	Niveau;		// Niveau d'importance de la trace.

public:
	int     Decal;          // decal


public:

	// Constructeur d'une trace.
	Trace(TYPE type, const string& nomClasse, const string& nomMethode, 
	      const string& message, void* pObjet= NULL,
	      NIVEAU niveau = IMPORTANT, int decal=0) : 
	      Type(type), NomClasse(nomClasse ), NomMethode(nomMethode),
              Message(message), Objet(pObjet) , Niveau(niveau), 
	      Decal(decal) {
	         INDENT.incIndent(Decal);   
	      }

	// Insertion d'une trace dans un flux
	virtual void InsererDansFlux(ostream& flux) const {
	     flux << ""		<< NomClasse;
	     flux << "::"	<< NomMethode;
	     flux << " : "	<< Message;
        }

	friend ostream& operator << (ostream& flux, const Trace& trace)	{
		trace.InsererDansFlux(flux);
		return flux;
	}
	
	virtual ~Trace() {
	   INDENT.decIndent(Decal);
	}
};







class Traceur  {

private:
	Traceur() : NiveauNormal(Trace::IMPORTANT),          // constructeur privé.
	            TypeNormal(Trace::METIER) 
		       {Screen = false; File = false;}	 
            
	static Traceur* Singleton;           // Unique Objet de la classe.

	Trace::NIVEAU	NiveauNormal;	     // Niveau maxi de trace à tracer
	Trace::TYPE	TypeNormal;          // Type de trace max à tracer
	ofstream	Flux;	             // Flux de sortie des trace.
	bool            Screen;              // screen Flux
	bool            File;                // file Flux	

public:

// Accesseurs 
	void Tracer(Trace* trace) {
	
	   if (trace->Niveau == Trace::IMPORTANT || trace->Type == Trace::ERROR) {
	        for (int i=0;i<INDENT.getIndent();i++) cout << " ";	
		if (Screen) cout << *trace << endl;
		if (File) Flux << *trace << endl;
		if (trace->Type == Trace::ERROR) exit(-1);
		return;
	   }

	   if (trace->Type <= TypeNormal && trace->Niveau <= NiveauNormal ) {
	        for (int i=0;i<INDENT.getIndent();i++) cout << " ";
		if (Screen) cout << (*trace) << endl;
		if (File) Flux << *trace << endl;
		return;
	   }
	};
	void Tracer(Trace& trace) {
		
	   if (trace.Niveau == Trace::IMPORTANT || trace.Type == Trace::ERROR) {
		for (int i=0;i<INDENT.getIndent();i++) cout << " ";	
		if (Screen) cout << trace << endl;
		if (File) Flux << trace << endl;
		if (trace.Type == Trace::ERROR) exit(-1);
		return;
	   }

	   if (trace.Type <= TypeNormal && trace.Niveau <= NiveauNormal ) {
		for (int i=0;i<INDENT.getIndent();i++) cout << " ";
		if (Screen) cout << trace << endl;
		if (File) Flux << trace << endl;
		return;
	   }
	};		
	void OuvrirFichierDeTraces(string nomFichier) {
	   if (Flux.is_open())
		Flux.close();
	   Flux.open(nomFichier.c_str());
	};	

	void FermerFichierDeTraces() {if (Flux.is_open()) Flux.close();};

	void SetNiveauTrace(Trace::NIVEAU niveau) 
		{NiveauNormal = niveau;};
	void SetTypeTrace(Trace::TYPE type) 
		{TypeNormal = type;};
	void SetScreenFlux (bool Flag)
	        {Screen=Flag;};
	void SetFileFlux (bool Flag, char* FileName)
	        {File=Flag; 
		if (File) OuvrirFichierDeTraces(FileName);};		
	void TraceMaximum()
		{NiveauNormal = Trace::DEBUG;
		 TypeNormal = Trace::TECHNIQUE;};		 
	void TraceMinimum() 
		{NiveauNormal = Trace::IMPORTANT;
		 TypeNormal = Trace::ERROR;};
	void TraceVerbose()
		{NiveauNormal = Trace::VERBOSE;
		 TypeNormal = Trace::METIER;};
	void TraceDebug()		 
		{NiveauNormal = Trace::DEBUG;
		 TypeNormal = Trace::MEMORY;};		 
		 
public :
        // Méthodes d'interface 
	void operator << (Trace* trace){	// Opérateur de flux.
		Tracer(trace);
	}
	void operator << (Trace& trace){	// Opérateur de flux.
		Tracer(trace);
	}	
        // Pattern Singleton  
	static Traceur& GetSingleton() {       // Accesseur statique
           if (Singleton != NULL) return *Singleton;
           Singleton = new Traceur();
           return *Singleton;
        };	        \
}; // Traceur



#endif
