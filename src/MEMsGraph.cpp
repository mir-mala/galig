#include "MEMsGraph.hpp"

MemsGraph::MemsGraph(const std::string& read_, const int& L_, const int& eps_, const int& exsN_, const bool& verbose_) : AnnNodesMap(AnnGraph),
                                             AnnEdgesMap(AnnGraph),
                                             NovNodesMap(NovGraph),
                                             NovEdgesMap(NovGraph) {
    read = read_;
    m = read.size();
    L = L_;
    eps = eps_;
    float t1 = eps*m*L;
    float t2 = eps*m;
    K0 = ceil(t1/100);
    K1 = ceil(t1/100) - L + 1;
    K2 = ceil(t2/100);
    //std::cout << m << " " << t2 << " " << K2 << std::endl;
    exsN = exsN_;
    verbose = verbose_;

    AnnStart = AnnGraph.addNode();
    AnnEnd = AnnGraph.addNode();
    NovStart = NovGraph.addNode();
    NovEnd = NovGraph.addNode();

    AnnNodesMap[AnnStart] = Mem(0,0,0);
    AnnNodesMap[AnnEnd] = Mem(-1,-1,-1);
    NovNodesMap[NovStart] = Mem(0,0,0);
    NovNodesMap[NovEnd] = Mem(-1,-1,-1);
}

std::pair<bool, int> MemsGraph::checkMEMs(const SplicingGraph& sg, const Mem& m1, const Mem& m2) {
    int id1 = sg.rank(m1.t-1);
    int id2 = sg.rank(m2.t-1);

    std::string exon1_text = sg.getExon(id1);
    std::string exon2_text = sg.getExon(id2);

    int err = -1;
    bool type = true;
    if(verbose) {
        std::cout << "Extending " << m1.toStr() << " with " << m2.toStr() << std::endl;
    }
    if(id1 == id2) { //m1 and m2 in the same exon
        //if(m2.p+m2.l>m1.p+m1.l && m1.t<m2.t && m2.t<m1.t+m1.l+K1 && m1.t+m1.l<m2.t+m2.l) {
        if(m2.p+m2.l>m1.p+m1.l && m1.t<m2.t && m1.t+m1.l<m2.t+m2.l) {
            if(verbose) {
                std::cout << "same exon" << std::endl;
            }
            int gapP = m2.p-m1.p-m1.l;
            int gap_E = m2.t-m1.t-m1.l;
            if(gapP>=0 && gap_E>=0) {
                if(gapP == 0) {
                    if(gap_E > K2) {
                        //Possible intron
                        if(verbose) {
                            std::cout << "Possible intron without overlap" << std::endl;
                        }
                        err = 0;
                        type = false;
                    } else {
                        //Errors
                        if(verbose) {
                            std::cout << "Nothing" << std::endl;
                        }
                        err = gap_E;
                        type = true;
                    }
                } else if(abs(gapP-gap_E) <= K2) {
                    //Possible SNV
                    if(verbose) {
                        std::cout << "Nothing" << std::endl;
                    }
                    std::string sub_P = read.substr(m1.p+m1.l-1, m2.p-m1.p-m1.l);
                    std::string sub_E = exon1_text.substr(m1.t+m1.l-sg.select(id1)-1-1, m2.t-m1.t-m1.l);
                    err = editDistance(sub_P, sub_E);
                    type = true;
                }
            } else if(gapP<=0 && gap_E<=0) {
                if(verbose) {
                    std::cout << "Nothing" << std::endl;
                }
                err = abs(gapP-gap_E);
                type = true;
            } else if(gapP<=0 && gap_E>K2) {
                //Possible intron
                if(verbose) {
                    std::cout << "Possible intron with overlap" << std::endl;
                }
                err = 0;
                type = false;
            } else {
                if(verbose) {
                    std::cout << "Nothing" << std::endl;
                }
                err = abs(gapP) + abs(gap_E);
                type = true;
            }
        }
    } else { //m1 and m2 in different exons
        if(sg.contain(id1, id2)) {
            if(verbose) {
                std::cout << "different exons" << std::endl;
            }
            if(m2.p+m2.l>m1.p+m1.l) {
                int gapP = m2.p-m1.p-m1.l;
                int gapE1 = sg.select(id1+1)+1-m1.t-m1.l;
                int gapE2 = m2.t-sg.select(id2)-1-1;
                if(gapP <= 0) {
                    err = 0; //abs(gapP);
                    if(verbose) {
                        std::cout << id1 << " " << id2 << " " << sg.isNew(id1, id2) << " " << gapE1 << " " << gapE2 << std::endl;
                    }
                    if(!sg.isNew(id1, id2) && gapE1 == 0 && gapE2 == 0) {
                        type = true;
                    }
                    else if(err <= K2) {
                        //Possible Competing
                        type = false;
                    } else
                        err = -1;
                } else {
                    if(gapE1 == 0 && gapE2 == 0) {
                        //Possible insertion (only if annotated edge)
                        if(!sg.isNew(id1,id2)) {
                            err = 0;
                            type = false;
                        }
                    } else {
                        if(abs(gapP-(gapE1+gapE2)) <= K2) {
                            //Possible SNV
                            if(verbose) {
                                std::cout << "SNV" << std::endl;
                            }
                            std::string subP = read.substr(m1.p+m1.l-1, gapP);
                            std::string subE1 = exon1_text.substr(m1.t+m1.l-sg.select(id1)-1-1, gapE1);
                            std::string subE2 = exon2_text.substr(0, gapE2);
                            std::string subE = subE1 + subE2;
                            err = editDistance(subP, subE);
                            if(!sg.isNew(id1, id2)) {
                                type = true;
                            } else {
                                type = false;
                            }
                        }
                    }
                }
            }
        } else {
            if(verbose) {
                std::cout << "no edge" << std::endl;
            }
        }
    }
    if(err > K2) {
        err = -1;
    }
    if(verbose) {
        std::cout << type << " " << err << std::endl;
    }
    return std::make_pair(type, err);
}

std::pair<bool, int> MemsGraph::validStart(const SplicingGraph& sg, const Mem& Mem) {
    if(Mem.p <= K0) {
        int err = K2+1;
        if(Mem.p == 1) {
            err = 0;
        } else {
            int id = sg.rank(Mem.t-1);
            std::string exon_text = sg.getExon(id);
            int l = Mem.p-1;
            std::string sub_P = read.substr(0,l);
            std::string sub_E;
            int exon_pref_len = Mem.t - sg.select(id)-1 - 1;
            if(exon_pref_len < l) {
                int shared_pref_len = l-exon_pref_len;
                std::string exon_pref = exon_text.substr(0, exon_pref_len);
                err = l;
                std::list<int> parents = sg.getParents(id);
                for(std::list<int>::iterator it=parents.begin(); it!=parents.end(); ++it) {
                    /**
                     * We look only at the father of the node,
                     * IF he is long enough, we get its suffix;
                     * ELSE we get all its label (without going further checking all its parents)
                     **/
                    int par = *it;
                    std::string par_text = sg.getExon(par);
                    if(sg.select(par+1)-shared_pref_len-sg.select(par)-1>=0) {
                        sub_E = par_text.substr(sg.select(par+1)-shared_pref_len-sg.select(par)-1, shared_pref_len) + exon_pref;
                    } else {
                        sub_E = par_text + exon_pref;
                    }
                    int curr_err = editDistance(sub_P, sub_E);
                    if(curr_err < err) {
                        err = curr_err;
                    }
                }
            } else {
                sub_E = exon_text.substr(Mem.t-l-sg.select(id)-1-1, l);
                err = editDistance(sub_P, sub_E);
            }
        }
        if(err <= K2) {
            return std::make_pair(true, err);
        }
    }
    return std::make_pair(false,K2+1);
}

std::pair<bool, int> MemsGraph::validEnd(const SplicingGraph& sg, const Mem& Mem) {
    if(Mem.p+Mem.l>=m-K0) {
        int err = K2+1;
        if(Mem.p+Mem.l == m+1) {
            err = 0;
        } else {
            int id = sg.rank(Mem.t-1);
            std::string exon_text = sg.getExon(id);
            int l = m-(Mem.p+Mem.l)+1;
            std::string sub_P = read.substr(Mem.p+Mem.l-1, l);
            std::string sub_E;
            int exon_suff_len = sg.select(id+1) - (Mem.t+Mem.l) + 1;
            if(exon_suff_len < l) {
                std::list<int> sons = sg.getSons(id);
                int shared_suff_len = l-exon_suff_len;
                std::string exon_suff;
                if(exon_suff_len == 0) {
                    exon_suff = "";
                } else {
                    exon_suff = exon_text.substr(Mem.t+Mem.l-sg.select(id)-1-1, exon_suff_len);
                }
                err = l;
                for(std::list<int>::iterator it=sons.begin(); it!=sons.end(); ++it) {
                    int son = *it;
                    std::string son_text = sg.getExon(son);
                    sub_E = exon_suff + son_text.substr(0, shared_suff_len);
                    int curr_err = editDistance(sub_P, sub_E);
                    if(curr_err < err) {
                        err = curr_err;
                    }
                }
            } else {
                sub_E = exon_text.substr(Mem.t+Mem.l-sg.select(id)-1-1, l);
                err = editDistance(sub_P, sub_E);
            }
        }
        if(err <= K2) {
            return std::make_pair(true, err);
        }
    }
    return std::make_pair(false,K2+1);
}

// mio codice --------------------------------------------------------------------------------------------------------------------------
struct comparatore_MEMs
{
    bool operator()( Mem const & a, Mem const & b ) const 
    {
        return a.l > b.l;
    }
};

void MemsGraph::build(const SplicingGraph& sg, std::list<Mem>& MEMs_) {  
    std::vector<std::forward_list<Mem> > MEMs (m+1);
    std::list<Mem>& mem_massimi=MEMs_;
    for(const Mem& m : MEMs_)
    {
        int p = m.p;
        if(MEMs[p].empty()) // controllo che un elemento ci sia prima di controllare il massimo
        {
            MEMs[p].push_front(m);
        }
        else
        {
            if(m.l>=(MEMs[p].front()).l) // tengo il massimo in alto
            {
                MEMs[p].push_front(m);
            }
            else
            {       
                MEMs[p].insert_after(MEMs[p].begin(), m);
            }
        }
    }
    mem_massimi.sort(comparatore_MEMs());
    bool percorso=false;
    while(!percorso && mem_massimi.size()!=0) // fino a quando non trovo un percorso
    {
        bool nodo_start=false;
        bool nodo_1_not_new;
        bool nodo_2_not_new;
        bool nodo_3_not_new;
        Mem& nodo_massimo=mem_massimi.front();    
        Node AnnNode1;
        Node NovNode1;
        if(nodo_massimo.isNew) 
        {
            nodo_1_not_new=false;
            AnnNode1 = AnnGraph.addNode();
            NovNode1 = NovGraph.addNode();
            nodo_massimo.setAnnNode(AnnNode1);
            nodo_massimo.setNovNode(NovNode1);
            AnnNodesMap[AnnNode1] = nodo_massimo;
            NovNodesMap[NovNode1] = nodo_massimo; 
            std::pair<bool, int> start_info = validStart(sg, nodo_massimo);
            if(start_info.first) 
            {
                // rimane comunque il problema di avere dei nodi start che in realtà hanno dei padri e quindi non sono ideali
                Arc arc = AnnGraph.addArc(AnnStart,AnnNode1);
                AnnEdgesMap[arc] = start_info.second;
                arc = NovGraph.addArc(NovStart,NovNode1);
                NovEdgesMap[arc] = start_info.second;
                nodo_start=true;
            }            
        }
        else 
        {
            nodo_1_not_new=true;
            AnnNode1 = nodo_massimo.AnnNode;
            NovNode1 = nodo_massimo.NovNode;
        }        
        int p2 = nodo_massimo.p+1;
        Mem figlio_massimo=mem_massimi.back();
        Mem figlio_massimo_ann=mem_massimi.back();
        Mem figlio_attuale=mem_massimi.back();
        bool figlio_massimo_ann_stato=false;   
        bool figlio_massimo_instanziato=false;
        bool figlio_stato=false; // se è true sono sicuro di avere un figlio POSSIBILE (err>=0)
        while(p2 <= m-L+1)
        {            
            if(!MEMs[p2].empty())
            {
                if(figlio_massimo_instanziato==false)
                {
                    figlio_massimo=MEMs[p2].front();
                    figlio_massimo_instanziato=true;
                }
                figlio_attuale=MEMs[p2].front();
                std::pair<bool, int> linkageInfo = checkMEMs(sg, nodo_massimo, figlio_attuale);
                bool flag = linkageInfo.first;
                int err = linkageInfo.second;
                if(err>=0) // se è vero allora questo è un mio possibile figlio
                {               
                    if(figlio_stato==false && figlio_attuale.l<=figlio_massimo.l) 
                    {
                        // se il primo MEMs[p2].front() ha err<0 ma ha il valore di l maggiore è un problema
                        figlio_massimo=figlio_attuale;
                        if(flag)
                        {
                            figlio_massimo_ann=figlio_attuale;
                            figlio_massimo_ann_stato=true;
                        }
                    }
                    else
                    {
                        if(figlio_attuale.l>figlio_massimo.l)
                        {
                            figlio_massimo=figlio_attuale;
                            if(flag) // può essere un figlio annotato
                            {
                                figlio_massimo_ann=figlio_attuale;
                                figlio_massimo_ann_stato=true;
                            }
                        }         
                        else
                        {
                            if(flag)
                            {
                                if(figlio_massimo_ann_stato==false) // per assicurarmi il figlio_massimo_ann
                                {
                                    figlio_massimo_ann=figlio_attuale;
                                    figlio_massimo_ann_stato=true;
                                }
                                else
                                {
                                    if(figlio_attuale.l>figlio_massimo_ann.l) // tengo anche il figlio maggiore annotato
                                    {
                                        figlio_massimo_ann=figlio_attuale;
                                    }
                                }
                            }
                        }    
                    }     
                    figlio_stato=true; // segno che ho un figlio             
                }
            }
            p2=p2+1;
        }
        Node AnnNode2;
        Node NovNode2;  
        Node NovNode3; // in caso il figlio massimo sia soltanto novel, ma ho trovato un figlio che comunque può essere annotato
        std::pair<bool, int> linkageInfo = checkMEMs(sg, nodo_massimo, figlio_massimo);
        bool flag = linkageInfo.first;
        int err = linkageInfo.second;
        if(figlio_stato)
        {
            if(figlio_massimo.isNew)
            {
                nodo_2_not_new=false;
                if(flag) // vuoldire che il figlio_massimo è uguale al figlio_massimo_ann
                {
                    AnnNode2 = AnnGraph.addNode();
                    NovNode2 = NovGraph.addNode();
                    figlio_massimo.setAnnNode(AnnNode2);
                    figlio_massimo.setNovNode(NovNode2);
                    AnnNodesMap[AnnNode2] = figlio_massimo;
                    NovNodesMap[NovNode2] = figlio_massimo;              
                }  
                else // aggiungo comunque un nodo annotato e due nodi novel per poter proseguire da quel nodo
                {
                    if(figlio_massimo_ann_stato)
                    {
                        if(figlio_massimo_ann.isNew)
                        {
                            nodo_3_not_new=false;
                            AnnNode2 = AnnGraph.addNode();
                            NovNode2 = NovGraph.addNode();
                            NovNode3 = NovGraph.addNode();
                            figlio_massimo_ann.setAnnNode(AnnNode2);
                            figlio_massimo_ann.setNovNode(NovNode2);
                            figlio_massimo.setNovNode(NovNode3);
                            AnnNodesMap[AnnNode2] = figlio_massimo_ann;
                            NovNodesMap[NovNode2] = figlio_massimo_ann;    
                            NovNodesMap[NovNode3] = figlio_massimo;  
                        }
                        else
                        {
                            nodo_3_not_new=true;
                            NovNode3 = NovGraph.addNode();
                            figlio_massimo.setNovNode(NovNode3); 
                            NovNodesMap[NovNode3] = figlio_massimo;  
                            AnnNode2 = figlio_massimo_ann.AnnNode;
                            NovNode2 = figlio_massimo_ann.NovNode;
                        }  
                    }        
                    else // non ho nessun figlio annotato, segno solo il massimo figlio novel
                    {                  
                        NovNode3 = NovGraph.addNode();
                        figlio_massimo.setNovNode(NovNode3);   
                        NovNodesMap[NovNode3] = figlio_massimo;  
                    }
                    
                }          
            } 
            else
            {
                nodo_2_not_new=true;
                if(flag)
                {
                    AnnNode2 = figlio_massimo.AnnNode;
                    NovNode2 = figlio_massimo.NovNode;
                }
                else
                {
                    if(figlio_massimo_ann_stato)
                    {
                        if(figlio_massimo_ann.isNew)
                        {
                            nodo_3_not_new=false;
                            AnnNode2 = AnnGraph.addNode();
                            NovNode2 = NovGraph.addNode();
                            figlio_massimo_ann.setAnnNode(AnnNode2);
                            figlio_massimo_ann.setNovNode(NovNode2);
                            AnnNodesMap[AnnNode2] = figlio_massimo_ann;
                            NovNodesMap[NovNode2] = figlio_massimo_ann;
                            NovNode3 = figlio_massimo.NovNode;    
                        }
                        else
                        {
                            nodo_3_not_new=true;
                            AnnNode2 = figlio_massimo_ann.AnnNode;
                            NovNode2 = figlio_massimo_ann.NovNode;
                            NovNode3 = figlio_massimo.NovNode;
                        }                        
                    }
                    else
                    {
                        NovNode3 = figlio_massimo.NovNode;
                    }
                }         
            }
            if(flag)
            {
                Arc arc = AnnGraph.addArc(AnnNode1,AnnNode2);
                AnnEdgesMap[arc] = err;
                arc = NovGraph.addArc(NovNode1,NovNode2);
                NovEdgesMap[arc] = err;
            }
            else 
            {
                if(figlio_massimo_ann_stato)
                {
                    Arc arc = AnnGraph.addArc(AnnNode1,AnnNode2);
                    AnnEdgesMap[arc] = err;
                    arc = NovGraph.addArc(NovNode1,NovNode2);
                    NovEdgesMap[arc] = err;
                    arc = NovGraph.addArc(NovNode1,NovNode3);
                    NovEdgesMap[arc] = err;
                }
                else
                {
                    Arc arc = NovGraph.addArc(NovNode1,NovNode3);
                    NovEdgesMap[arc] = err;
                }
            }
        }
        else // se ha figli non ha senso che sia un nodo di end
        {
            std::pair<bool, int> end_info = validEnd(sg, nodo_massimo); // creo la lista dei nodi end
            if(end_info.first)
            {
                Arc arc = AnnGraph.addArc(AnnNode1,AnnEnd);
                AnnEdgesMap[arc] = end_info.second;
                arc = NovGraph.addArc(NovNode1,NovEnd);
                NovEdgesMap[arc] = end_info.second;
            }
        }
        Node nodo_partenza;
        bool controlla=false;
        if(figlio_stato) // altrimenti non ho fatto collegamenti che possono aver creato un path
        {
            if(flag) // caso in cui il figlio massimo è anche annotato
            {
                if(nodo_1_not_new==true && nodo_2_not_new==true) // cerco se esiste un path
                {                 
                    controlla=true;
                    nodo_partenza=AnnStart;
                }
                // if(nodo_1_not_new==true && nodo_2_not_new==false) non ha senso controllarlo perchè dovrei controllare
                // in modo ricorsivo se il secondo nodo ha figli fino ad arrivare ad un nodo di end,
                // il secondo nodo non mi garantisce di aver costruito un percorso
                if(nodo_1_not_new==false && nodo_2_not_new==true && nodo_start==true)
                {
                    // cerco se esiste il path partendo da questo particolare nodo di start
                    controlla=true;
                    nodo_partenza=AnnNode1;                    
                }
            }
            else // controllo di avere un figlio massimo annotato anche se non è il vero figlio massimo
            {                
                if(figlio_massimo_ann_stato)
                {
                    if(nodo_1_not_new==true && nodo_3_not_new==true)
                    {
                        controlla=true;
                        nodo_partenza=AnnStart;                         
                    }
                    if(nodo_1_not_new==false && nodo_3_not_new==true && nodo_start==true)
                    {
                        controlla=true;
                        nodo_partenza=AnnNode1;  
                    }
                }  
            }
        }
        if(controlla)
        {
            lemon::Dijkstra<Graph, lemon::ListDigraph::ArcMap<int> >
                        ::SetStandardHeap<FibH>
                        ::SetHeap<FibH,FibM>
                        ::Create AnnDijkstra (AnnGraph, AnnEdgesMap);
            FibM AnnHCR (AnnGraph);
            FibH AnnHeap (AnnHCR);
            AnnDijkstra.heap(AnnHeap, AnnHCR);
            AnnDijkstra.run(nodo_partenza,AnnEnd);
            if(AnnDijkstra.reached(AnnEnd))
            {
                percorso=true;;    
            }
        }
        mem_massimi.pop_front();
    }    
    save("Graph.dot");   
}
// fine mio codice --------------------------------------------------------------------------------------------------------------------------

std::list<std::pair<int, std::list<Mem> > > MemsGraph::visit(const SplicingGraph& sg) {
    std::list<std::pair<int, std::list<Mem> > > paths;
    std::list<Mem> AnnPath1;
    std::list<Mem> AnnPath2;
    std::list<Mem> NovPath;
    bool FoundAnnotated = false;
    int AnnW1 = K2+1;
    int AnnW2 = K2+1;
    int NovW = K2+1;

    //Visiting Annotated Graph
    lemon::Dijkstra<Graph, lemon::ListDigraph::ArcMap<int> >
        ::SetStandardHeap<FibH>
        ::SetHeap<FibH,FibM>
        ::Create AnnDijkstra (AnnGraph, AnnEdgesMap);
    FibM AnnHCR (AnnGraph);
    FibH AnnHeap (AnnHCR);
    AnnDijkstra.heap(AnnHeap, AnnHCR);
    AnnDijkstra.run(AnnStart,AnnEnd);
    if(AnnDijkstra.reached(AnnEnd)) {
        AnnW1 = AnnDijkstra.dist(AnnEnd);
        if(AnnW1 <= K2) {
            FoundAnnotated = true;
            Path p = AnnDijkstra.path(AnnEnd);
            bool first = true;
            for(Path::ArcIt it(p); it != lemon::INVALID; ++it) {
                if(first) {
                    AnnGraph.erase(it);
                    first = false;
                }
                Arc e = it;
                Node target = AnnGraph.target(e);
                Mem m = AnnNodesMap[target];
                if(NovGraph.id(target) != NovGraph.id(AnnEnd)) {
                    AnnPath1.push_back(m);
                }
            }
            paths.push_back(std::make_pair(AnnW1, AnnPath1));
        }
    }

    AnnDijkstra.run(AnnStart,AnnEnd);
    if(AnnDijkstra.reached(AnnEnd)) {
        AnnW2 = AnnDijkstra.dist(AnnEnd);
        if(AnnW2 <= K2) {
            Path p = AnnDijkstra.path(AnnEnd);
            for(Path::ArcIt it(p); it != lemon::INVALID; ++it) {
                Arc e = it;
                Node target = AnnGraph.target(e);
                Mem m = AnnNodesMap[target];
                if(NovGraph.id(target) != NovGraph.id(AnnEnd)) {
                    AnnPath2.push_back(m);
                }
            }
            paths.push_back(std::make_pair(AnnW2, AnnPath2));
        }
    }

    //Visiting Novel Graph
    lemon::Dijkstra<Graph, lemon::ListDigraph::ArcMap<int> >
        ::SetStandardHeap<FibH>
        ::SetHeap<FibH,FibM>
        ::Create NovDijkstra (NovGraph, NovEdgesMap);
    FibM NovHCR (NovGraph);
    FibH NovHeap (NovHCR);
    NovDijkstra.heap(NovHeap, NovHCR);

    NovDijkstra.run(NovStart,NovEnd);
    if(NovDijkstra.reached(NovEnd)) {
        NovW = NovDijkstra.dist(NovEnd);
        if(!FoundAnnotated) {
            if(NovW <= K2) {
                Path p = NovDijkstra.path(NovEnd);
                for(Path::ArcIt it(p); it != lemon::INVALID; ++it) {
                    Arc e = it;
                    Node target = NovGraph.target(e);
                    Mem m = NovNodesMap[target];
                    if(NovGraph.id(target) != NovGraph.id(NovEnd)) {
                        NovPath.push_back(m);
                    }
                }
            }
            if(paths.size() > 1)
                paths.pop_back();
            paths.push_front(std::make_pair(NovW, NovPath));
        }
    }
    return paths;
}

void MemsGraph::save(const std::string& s) {
    std::ofstream myfile;

    myfile.open("Ann"+s);

    std::string dot = "digraph G {\n graph [splines=true overlap=false]\n node  [shape=ellipse, width=0.3, height=0.3]\n";
    for(NodeIt n (AnnGraph); n != lemon::INVALID; ++n) {
        dot += " " + std::to_string(AnnGraph.id(n)) + " [label=\"" + AnnNodesMap[n].toStr() + "\"];\n";
    }
    for(ArcIt a (AnnGraph); a != lemon::INVALID; ++a) {
        dot += " " + std::to_string(AnnGraph.id(AnnGraph.source(a))) + " -> " + std::to_string(AnnGraph.id(AnnGraph.target(a))) + "[label=\"" + std::to_string(AnnEdgesMap[a]) + "\"];\n";
    }
    dot += "}";

    myfile << dot;
    myfile.close();
    
    myfile.open("Nov"+s);

    dot = "digraph G {\n graph [splines=true overlap=false]\n node  [shape=ellipse, width=0.3, height=0.3]\n";
    for (NodeIt n (NovGraph); n != lemon::INVALID; ++n) {
        dot += " " + std::to_string(NovGraph.id(n)) + " [label=\"" + NovNodesMap[n].toStr() + "\"];\n";
    }
    for(ArcIt a (NovGraph); a != lemon::INVALID; ++a) {
        dot += " " + std::to_string(NovGraph.id(NovGraph.source(a))) + " -> " + std::to_string(NovGraph.id(NovGraph.target(a))) + "[label=\"" + std::to_string(NovEdgesMap[a]) + "\"];\n";
    }
    dot += "}";

    myfile << dot;
    myfile.close();
}
