#include "MEMsGraph.hpp"

MemsGraph::MemsGraph(const std::string &read_, const int &L_, const int &eps_, const int &exsN_, const bool &verbose_) : AnnNodesMap(AnnGraph),
                                                                                                                         AnnEdgesMap(AnnGraph),
                                                                                                                         NovNodesMap(NovGraph),
                                                                                                                         NovEdgesMap(NovGraph)
{
    read = read_;
    m = read.size();
    L = L_;
    eps = eps_;
    float t1 = eps * m * L;
    float t2 = eps * m;
    K0 = ceil(t1 / 100);
    K1 = ceil(t1 / 100) - L + 1;
    K2 = ceil(t2 / 100);
    //std::cout << m << " " << t2 << " " << K2 << std::endl;
    exsN = exsN_;
    verbose = verbose_;

    AnnStart = AnnGraph.addNode();
    AnnEnd = AnnGraph.addNode();
    NovStart = NovGraph.addNode();
    NovEnd = NovGraph.addNode();

    AnnNodesMap[AnnStart] = Mem(0, 0, 0);
    AnnNodesMap[AnnEnd] = Mem(-1, -1, -1);
    NovNodesMap[NovStart] = Mem(0, 0, 0);
    NovNodesMap[NovEnd] = Mem(-1, -1, -1);
}

std::pair<bool, int> MemsGraph::checkMEMs(const SplicingGraph &sg, const Mem &m1, const Mem &m2)
{
    int id1 = sg.rank(m1.t - 1);
    int id2 = sg.rank(m2.t - 1);

    std::string exon1_text = sg.getExon(id1);
    std::string exon2_text = sg.getExon(id2);

    int err = -1;
    bool type = true;
    if (verbose)
    {
        std::cout << "Extending " << m1.toStr() << " with " << m2.toStr() << std::endl;
    }
    if (id1 == id2)
    { //m1 and m2 in the same exon
        //if(m2.p+m2.l>m1.p+m1.l && m1.t<m2.t && m2.t<m1.t+m1.l+K1 && m1.t+m1.l<m2.t+m2.l) {
        if (m2.p + m2.l > m1.p + m1.l && m1.t < m2.t && m1.t + m1.l < m2.t + m2.l)
        {
            if (verbose)
            {
                std::cout << "same exon" << std::endl;
            }
            int gapP = m2.p - m1.p - m1.l;
            int gap_E = m2.t - m1.t - m1.l;
            if (gapP >= 0 && gap_E >= 0)
            {
                if (gapP == 0)
                {
                    if (gap_E > K2)
                    {
                        //Possible intron
                        if (verbose)
                        {
                            std::cout << "Possible intron without overlap" << std::endl;
                        }
                        err = 0;
                        type = false;
                    }
                    else
                    {
                        //Errors
                        if (verbose)
                        {
                            std::cout << "Nothing" << std::endl;
                        }
                        err = gap_E;
                        type = true;
                    }
                }
                else if (abs(gapP - gap_E) <= K2)
                {
                    //Possible SNV
                    if (verbose)
                    {
                        std::cout << "Nothing" << std::endl;
                    }
                    std::string sub_P = read.substr(m1.p + m1.l - 1, m2.p - m1.p - m1.l);
                    std::string sub_E = exon1_text.substr(m1.t + m1.l - sg.select(id1) - 1 - 1, m2.t - m1.t - m1.l);
                    err = editDistance(sub_P, sub_E);
                    type = true;
                }
            }
            else if (gapP <= 0 && gap_E <= 0)
            {
                if (verbose)
                {
                    std::cout << "Nothing" << std::endl;
                }
                err = abs(gapP - gap_E);
                type = true;
            }
            else if (gapP <= 0 && gap_E > K2)
            {
                //Possible intron
                if (verbose)
                {
                    std::cout << "Possible intron with overlap" << std::endl;
                }
                err = 0;
                type = false;
            }
            else
            {
                if (verbose)
                {
                    std::cout << "Nothing" << std::endl;
                }
                err = abs(gapP) + abs(gap_E);
                type = true;
            }
        }
    }
    else
    { //m1 and m2 in different exons
        if (sg.contain(id1, id2))
        {
            if (verbose)
            {
                std::cout << "different exons" << std::endl;
            }
            if (m2.p + m2.l > m1.p + m1.l)
            {
                int gapP = m2.p - m1.p - m1.l;
                int gapE1 = sg.select(id1 + 1) + 1 - m1.t - m1.l;
                int gapE2 = m2.t - sg.select(id2) - 1 - 1;
                if (gapP <= 0)
                {
                    err = 0; //abs(gapP);
                    if (verbose)
                    {
                        std::cout << id1 << " " << id2 << " " << sg.isNew(id1, id2) << " " << gapE1 << " " << gapE2 << std::endl;
                    }
                    if (!sg.isNew(id1, id2) && gapE1 == 0 && gapE2 == 0)
                    {
                        type = true;
                    }
                    else if (err <= K2)
                    {
                        //Possible Competing
                        type = false;
                    }
                    else
                        err = -1;
                }
                else
                {
                    if (gapE1 == 0 && gapE2 == 0)
                    {
                        //Possible insertion (only if annotated edge)
                        if (!sg.isNew(id1, id2))
                        {
                            err = 0;
                            type = false;
                        }
                    }
                    else
                    {
                        if (abs(gapP - (gapE1 + gapE2)) <= K2)
                        {
                            //Possible SNV
                            if (verbose)
                            {
                                std::cout << "SNV" << std::endl;
                            }
                            std::string subP = read.substr(m1.p + m1.l - 1, gapP);
                            std::string subE1 = exon1_text.substr(m1.t + m1.l - sg.select(id1) - 1 - 1, gapE1);
                            std::string subE2 = exon2_text.substr(0, gapE2);
                            std::string subE = subE1 + subE2;
                            err = editDistance(subP, subE);
                            if (!sg.isNew(id1, id2))
                            {
                                type = true;
                            }
                            else
                            {
                                type = false;
                            }
                        }
                    }
                }
            }
        }
        else
        {
            if (verbose)
            {
                std::cout << "no edge" << std::endl;
            }
        }
    }
    if (err > K2)
    {
        err = -1;
    }
    if (verbose)
    {
        std::cout << type << " " << err << std::endl;
    }
    return std::make_pair(type, err);
}

std::pair<bool, int> MemsGraph::validStart(const SplicingGraph &sg, const Mem &Mem)
{
    if (Mem.p <= K0)
    {
        int err = K2 + 1;
        if (Mem.p == 1)
        {
            err = 0;
        }
        else
        {
            int id = sg.rank(Mem.t - 1);
            std::string exon_text = sg.getExon(id);
            int l = Mem.p - 1;
            std::string sub_P = read.substr(0, l);
            std::string sub_E;
            int exon_pref_len = Mem.t - sg.select(id) - 1 - 1;
            if (exon_pref_len < l)
            {
                int shared_pref_len = l - exon_pref_len;
                std::string exon_pref = exon_text.substr(0, exon_pref_len);
                err = l;
                std::list<int> parents = sg.getParents(id);
                for (std::list<int>::iterator it = parents.begin(); it != parents.end(); ++it)
                {
                    /**
                     * We look only at the father of the node,
                     * IF he is long enough, we get its suffix;
                     * ELSE we get all its label (without going further checking all its parents)
                     **/
                    int par = *it;
                    std::string par_text = sg.getExon(par);
                    if (sg.select(par + 1) - shared_pref_len - sg.select(par) - 1 >= 0)
                    {
                        sub_E = par_text.substr(sg.select(par + 1) - shared_pref_len - sg.select(par) - 1, shared_pref_len) + exon_pref;
                    }
                    else
                    {
                        sub_E = par_text + exon_pref;
                    }
                    int curr_err = editDistance(sub_P, sub_E);
                    if (curr_err < err)
                    {
                        err = curr_err;
                    }
                }
            }
            else
            {
                sub_E = exon_text.substr(Mem.t - l - sg.select(id) - 1 - 1, l);
                err = editDistance(sub_P, sub_E);
            }
        }
        if (err <= K2)
        {
            return std::make_pair(true, err);
        }
    }
    return std::make_pair(false, K2 + 1);
}

std::pair<bool, int> MemsGraph::validEnd(const SplicingGraph &sg, const Mem &Mem)
{
    if (Mem.p + Mem.l >= m - K0)
    {
        int err = K2 + 1;
        if (Mem.p + Mem.l == m + 1)
        {
            err = 0;
        }
        else
        {
            int id = sg.rank(Mem.t - 1);
            std::string exon_text = sg.getExon(id);
            int l = m - (Mem.p + Mem.l) + 1;
            std::string sub_P = read.substr(Mem.p + Mem.l - 1, l);
            std::string sub_E;
            int exon_suff_len = sg.select(id + 1) - (Mem.t + Mem.l) + 1;
            if (exon_suff_len < l)
            {
                std::list<int> sons = sg.getSons(id);
                int shared_suff_len = l - exon_suff_len;
                std::string exon_suff;
                if (exon_suff_len == 0)
                {
                    exon_suff = "";
                }
                else
                {
                    exon_suff = exon_text.substr(Mem.t + Mem.l - sg.select(id) - 1 - 1, exon_suff_len);
                }
                err = l;
                for (std::list<int>::iterator it = sons.begin(); it != sons.end(); ++it)
                {
                    int son = *it;
                    std::string son_text = sg.getExon(son);
                    sub_E = exon_suff + son_text.substr(0, shared_suff_len);
                    int curr_err = editDistance(sub_P, sub_E);
                    if (curr_err < err)
                    {
                        err = curr_err;
                    }
                }
            }
            else
            {
                sub_E = exon_text.substr(Mem.t + Mem.l - sg.select(id) - 1 - 1, l);
                err = editDistance(sub_P, sub_E);
            }
        }
        if (err <= K2)
        {
            return std::make_pair(true, err);
        }
    }
    return std::make_pair(false, K2 + 1);
}

// mio codice --------------------------------------------------------------------------------------------------------------------------
struct comparatore_MEMs
{
    bool operator()(Mem const &a, Mem const &b) const
    {
        return a.l > b.l;
    }
};

//m1 è il padre o figlio massimo, m2 è il nodo massimo, m3 è il padre o figlio attuale, m4 è il padre o il filgio massimo_annotato
std::tuple<bool, bool, Mem, Mem> MemsGraph::cerca_nodo(int p2, int valore_confronto, const SplicingGraph &sg, Mem &m1, Mem &m2, Mem &m3, Mem &m4, bool ann_stato, bool stato,
                                                        std::vector<std::forward_list<Mem>> &MEMs)
{
    bool pf = false; //true sto trattando il padre, false sto trattando il figlio
    if (p2 == 0)
    {
        pf = true;
    }
    bool instanziato = false;
    while (p2 <= valore_confronto) //considero un overlap
    {
        if (!MEMs[p2].empty())
        {
            if (instanziato == false)
            {
                m1 = MEMs[p2].front();
                instanziato = true;
            }
            m3 = MEMs[p2].front();
            std::pair<bool, int> linkageInfo;
            if (pf == true)
            {
                linkageInfo = checkMEMs(sg, m3, m2);
            }
            else
            {
                linkageInfo = checkMEMs(sg, m2, m3);
            }
            bool flag = linkageInfo.first;
            int err = linkageInfo.second;
            if (err >= 0) // se è vero allora questo è un mio possibile (padre o figlio) 
            {
                if (stato == false)
                {
                    stato = true; // segno che ho un (padre o figlio)
                    // se il primo MEMs[p2].front() ha err<0 ma ha il valore di l maggiore è un problema
                    m1 = m3;
                    if (flag) // può essere un (padre o figlio) annotato
                    {
                        m4 = m3;
                        ann_stato = true;
                    }
                }
                else
                {
                    if (m3.l > m1.l)
                    {
                        m1 = m3;
                        if (flag)
                        {
                            m4 = m3;
                            ann_stato = true;
                        }
                    }
                    else
                    {
                        if (flag)
                        {
                            if (ann_stato == false) // per assicurarmi il (padre o figlio)_massimo_annotato
                            {
                                m4 = m3;
                                ann_stato = true;
                            }
                            else
                            {
                                if (m3.l > m4.l) // tengo aggiornato anche il (padre o figlio)_masssimo_annotato
                                {
                                    m4 = m3;
                                }
                            }
                        }
                    }
                }
            }
        }
        p2 = p2 + 1;
    }
    return std::make_tuple(ann_stato, stato, m1, m4);
}

// m1 è il padre o figlio massimo, m2 è il nodo massimo, m3 è il padre o figlio massimo_annotato, fp se true stiamo trattando il figlio, il padre altrimenti
std::tuple<bool, bool, Node, Node, Node> MemsGraph::collegamento_nodo(const SplicingGraph &sg, Mem &m1, Mem &m2, bool ann_stato, Mem &m3, Node AnnNode1, Node NovNode1, bool fp)
{
    Node AnnNode2;
    Node NovNode2;
    Node NovNode3; // in caso il (padre o figlio) massimo sia soltanto novel, ma ho trovato un (padre o figlio) che comunque può essere annotato
    std::pair<bool, int> linkageInfo;
    if (fp == true)
    {
        linkageInfo = checkMEMs(sg, m2, m1);
    }
    else
    {
        linkageInfo = checkMEMs(sg, m1, m2);
    }
    bool flag = linkageInfo.first;
    int err = linkageInfo.second;
    bool isnew = false;
    //creo il nodo
    if (m1.isNew)
    {
        isnew = true;
        if (flag) // vuoldire che il (padre o figlio)_massimo è uguale al (padre o figlio)_massimo_annotato
        {
            AnnNode2 = AnnGraph.addNode();
            NovNode2 = NovGraph.addNode();
            m1.setAnnNode(AnnNode2);
            m1.setNovNode(NovNode2);
            AnnNodesMap[AnnNode2] = m1;
            NovNodesMap[NovNode2] = m1;
        }
        else // aggiungo comunque un nodo annotato e due nodi novel per poter proseguire dal nodo
        {
            if (ann_stato)
            {
                if (m3.isNew)
                {
                    AnnNode2 = AnnGraph.addNode();
                    NovNode2 = NovGraph.addNode();
                    NovNode3 = NovGraph.addNode();
                    m3.setAnnNode(AnnNode2);
                    m3.setNovNode(NovNode2);
                    m1.setNovNode(NovNode3);
                    AnnNodesMap[AnnNode2] = m3;
                    NovNodesMap[NovNode2] = m3;
                    NovNodesMap[NovNode3] = m1;
                }
                else
                {
                    NovNode3 = NovGraph.addNode();
                    m1.setNovNode(NovNode3);
                    NovNodesMap[NovNode3] = m1;
                    AnnNode2 = m3.AnnNode;
                    NovNode2 = m3.NovNode;
                }
            }
            else // non ho nessun (padre o figlio) annotato, segno solo il massimo (padre o figlio) novel
            {
                NovNode3 = NovGraph.addNode();
                m1.setNovNode(NovNode3);
                NovNodesMap[NovNode3] = m1;
            }
        }
    }
    else
    {
        if (flag)
        {
            AnnNode2 = m1.AnnNode;
            NovNode2 = m1.NovNode;
        }
        else
        {
            if (ann_stato)
            {
                if (m3.isNew)
                {
                    AnnNode2 = AnnGraph.addNode();
                    NovNode2 = NovGraph.addNode();
                    m3.setAnnNode(AnnNode2);
                    m3.setNovNode(NovNode2);
                    AnnNodesMap[AnnNode2] = m3;
                    NovNodesMap[NovNode2] = m3;
                    NovNode3 = m1.NovNode;
                }
                else
                {
                    AnnNode2 = m3.AnnNode;
                    NovNode2 = m3.NovNode;
                    NovNode3 = m1.NovNode;
                }
            }
            else
            {
                NovNode3 = m1.NovNode;
            }
        }
    }
    //collego il nodo
    if (flag)
    {
        if (fp == true)
        {
            Arc arc = AnnGraph.addArc(AnnNode1, AnnNode2);
            AnnEdgesMap[arc] = err;
            arc = NovGraph.addArc(NovNode1, NovNode2);
            NovEdgesMap[arc] = err;
        }
        else
        {
            Arc arc = AnnGraph.addArc(AnnNode2, AnnNode1);
            AnnEdgesMap[arc] = err;
            arc = NovGraph.addArc(NovNode2, NovNode1);
            NovEdgesMap[arc] = err;
        }
    }
    else
    {
        if (ann_stato)
        {
            if (fp == true)
            {
                std::pair<bool, int> linkageInfo_ann = checkMEMs(sg, m2, m3);
                Arc arc = AnnGraph.addArc(AnnNode1, AnnNode2);
                AnnEdgesMap[arc] = linkageInfo_ann.second; // devo dargli l'errore corretto
                arc = NovGraph.addArc(NovNode1, NovNode2);
                NovEdgesMap[arc] = linkageInfo_ann.second;
                arc = NovGraph.addArc(NovNode1, NovNode3);
                NovEdgesMap[arc] = err;
            }
            else
            {
                std::pair<bool, int> linkageInfo_ann = checkMEMs(sg, m3, m2);
                Arc arc = AnnGraph.addArc(AnnNode2, AnnNode1);
                AnnEdgesMap[arc] = linkageInfo_ann.second; 
                arc = NovGraph.addArc(NovNode2, NovNode1);
                NovEdgesMap[arc] = linkageInfo_ann.second;
                arc = NovGraph.addArc(NovNode3, NovNode1);
                NovEdgesMap[arc] = err;
            }
        }
        else
        {
            if (fp == true)
            {
                Arc arc = NovGraph.addArc(NovNode1, NovNode3);
                NovEdgesMap[arc] = err;
            }
            else
            {
                Arc arc = NovGraph.addArc(NovNode3, NovNode1);
                NovEdgesMap[arc] = err;
            }
        }
    }
    return std::make_tuple(isnew, flag, AnnNode2, NovNode2, NovNode3);
}

void MemsGraph::build(const SplicingGraph &sg, std::list<Mem> &MEMs_)
{
    std::vector<std::forward_list<Mem>> MEMs(m + 1);
    std::list<Mem> &mem_massimi = MEMs_;
    for (const Mem &m : MEMs_)
    {
        int p = m.p;
        if (MEMs[p].empty()) // controllo che un elemento ci sia prima di controllare il massimo
        {
            MEMs[p].push_front(m);
        }
        else
        {
            if (m.l >= (MEMs[p].front()).l) // tengo il massimo in alto
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
    bool percorso = false;
    while (!percorso && mem_massimi.size() != 0) // fino a quando non trovo un percorso o visito tutti i MEM
    {
        Mem &nodo_massimo = mem_massimi.front();
        Node AnnNode1;
        Node NovNode1;
        bool nodo_massimo_isnew = nodo_massimo.isNew;
        bool padre_isnew = false;
        bool figlio_isnew = false;
        bool valid_start = false;
        bool valid_end = false;
        if (nodo_massimo_isnew)
        {
            AnnNode1 = AnnGraph.addNode();
            NovNode1 = NovGraph.addNode();
            nodo_massimo.setAnnNode(AnnNode1);
            nodo_massimo.setNovNode(NovNode1);
            AnnNodesMap[AnnNode1] = nodo_massimo;
            NovNodesMap[NovNode1] = nodo_massimo;
        }
        else
        {
            AnnNode1 = nodo_massimo.AnnNode;
            NovNode1 = nodo_massimo.NovNode;
        }
        // cerco il padre
        int p2 = 0;
        Mem padre_massimo = mem_massimi.back();
        Mem padre_massimo_ann = mem_massimi.back();
        Mem padre_attuale = mem_massimi.back();
        bool padre_massimo_ann_stato = false;
        bool padre_stato = false;                  // se è true sono sicuro di avere un padre POSSIBILE (err>=0)
        int valore_confronto = nodo_massimo.p - 1; // considero un overlap altrimenti dovrei fare "nodo_massimo.p-p2+1>=L"
        std::tie(padre_massimo_ann_stato, padre_stato, padre_massimo, padre_massimo_ann) = cerca_nodo(p2, valore_confronto, sg, padre_massimo, nodo_massimo, padre_attuale,
                                                                                                        padre_massimo_ann, padre_massimo_ann_stato, padre_stato, MEMs);
        // collego il nodo padre
        if (padre_stato)
        {
            bool flag;
            Node AnnNode2;
            Node NovNode2;
            Node NovNode3;
            std::tie(padre_isnew, flag, AnnNode2, NovNode2, NovNode3) = collegamento_nodo(sg, padre_massimo, nodo_massimo, padre_massimo_ann_stato, padre_massimo_ann, AnnNode1, NovNode1,
                                                                                            false);
            // controllo che il padre possa essere un nodo di start
            std::pair<bool, int> start_info_dad = validStart(sg, padre_massimo);
            if (flag)
            {
                if (start_info_dad.first)
                {
                    Arc arc = AnnGraph.addArc(AnnStart, AnnNode2);
                    AnnEdgesMap[arc] = start_info_dad.second;
                    arc = NovGraph.addArc(NovStart, NovNode2);
                    NovEdgesMap[arc] = start_info_dad.second;
                    valid_start = true;
                }
            }
            else
            {
                if (start_info_dad.first)
                {
                    Arc arc = NovGraph.addArc(NovStart, NovNode3);
                    NovEdgesMap[arc] = start_info_dad.second;
                    valid_start = true;
                }
                if (padre_massimo_ann_stato)
                {
                    std::pair<bool, int> start_info_dad_ann = validStart(sg, padre_massimo_ann);
                    if (start_info_dad_ann.first)
                    {
                        Arc arc = AnnGraph.addArc(AnnStart, AnnNode2);
                        AnnEdgesMap[arc] = start_info_dad_ann.second;
                        arc = NovGraph.addArc(NovStart, NovNode2);
                        NovEdgesMap[arc] = start_info_dad_ann.second;
                        valid_start = true;
                    }
                }
            }
        }
        else // controllo se il nodo massimo può essere un buon nodo di start
        {
            std::pair<bool, int> start_info_iniziale = validStart(sg, nodo_massimo);
            if (start_info_iniziale.first && nodo_massimo_isnew)
            {
                // è un buon nodo di start
                Arc arc = AnnGraph.addArc(AnnStart, AnnNode1);
                AnnEdgesMap[arc] = start_info_iniziale.second;
                arc = NovGraph.addArc(NovStart, NovNode1);
                NovEdgesMap[arc] = start_info_iniziale.second;
                valid_start = true;
            }
        }
        // cerco il figlio
        p2 = nodo_massimo.p + nodo_massimo.l - L + 1;
        ;
        Mem figlio_massimo = mem_massimi.back();
        Mem figlio_massimo_ann = mem_massimi.back();
        Mem figlio_attuale = mem_massimi.back();
        bool figlio_massimo_ann_stato = false;
        bool figlio_stato = false; // se è true sono sicuro di avere un figlio POSSIBILE (err>=0)
        valore_confronto = m - L + 1;
        std::tie(figlio_massimo_ann_stato, figlio_stato, figlio_massimo, figlio_massimo_ann) = cerca_nodo(p2, valore_confronto, sg, figlio_massimo, nodo_massimo, figlio_attuale,
                                                                                                            figlio_massimo_ann, figlio_massimo_ann_stato, figlio_stato, MEMs);
        //collego il figlio
        if (figlio_stato)
        {
            bool flag;
            Node AnnNode2;
            Node NovNode2;
            Node NovNode3;
            std::tie(figlio_isnew, flag, AnnNode2, NovNode2, NovNode3) = collegamento_nodo(sg, figlio_massimo, nodo_massimo, figlio_massimo_ann_stato, figlio_massimo_ann, AnnNode1,
                                                                                            NovNode1, true);
            // controllo che il figlio possa essere un nodo di end
            std::pair<bool, int> end_info_son = validEnd(sg, figlio_massimo);
            if (flag)
            {
                if (end_info_son.first)
                {
                    Arc arc = AnnGraph.addArc(AnnNode2, AnnEnd);
                    AnnEdgesMap[arc] = end_info_son.second;
                    arc = NovGraph.addArc(NovNode2, NovEnd);
                    NovEdgesMap[arc] = end_info_son.second;
                    valid_end = true;
                }
            }
            else
            {
                if (end_info_son.first)
                {
                    Arc arc = NovGraph.addArc(NovNode3, NovEnd);
                    NovEdgesMap[arc] = end_info_son.second;
                    valid_end = true;
                }
                if (figlio_massimo_ann_stato)
                {
                    std::pair<bool, int> end_info_son_ann = validEnd(sg, figlio_massimo_ann);
                    if (end_info_son_ann.first)
                    {
                        Arc arc = AnnGraph.addArc(AnnNode2, AnnEnd);
                        AnnEdgesMap[arc] = end_info_son_ann.second;
                        arc = NovGraph.addArc(NovNode2, NovEnd);
                        NovEdgesMap[arc] = end_info_son_ann.second;
                        valid_end = true;
                    }
                }
            }
        }
        else // controllo se il nodo massimo può essere un buon nodo di end
        {
            std::pair<bool, int> end_info_iniziale = validEnd(sg, nodo_massimo);
            if (end_info_iniziale.first && nodo_massimo_isnew)
            {
                // può essere un buon nodo di end
                Arc arc = AnnGraph.addArc(AnnNode1, AnnEnd);
                AnnEdgesMap[arc] = end_info_iniziale.second;
                arc = NovGraph.addArc(NovNode1, NovEnd);
                NovEdgesMap[arc] = end_info_iniziale.second;
                valid_end = true;
            }
        }
        if (valid_start && valid_end) // sono certo di aver trovato un percorso senza doverlo cercare
        {
            percorso = true;
        }
        bool controllo = false;
        if (percorso == false) // controllo la possibilità di aver trovato un percorso
        {
            if ((valid_start && !figlio_isnew) || (valid_end && !padre_isnew) || (!padre_isnew && !figlio_isnew))
            {
                controllo = true;
            }
        }
        if (controllo)
        {
            lemon::Dijkstra<Graph, lemon::ListDigraph::ArcMap<int>>::SetStandardHeap<FibH>::SetHeap<FibH, FibM>::Create AnnDijkstra(AnnGraph, AnnEdgesMap);
            FibM AnnHCR(AnnGraph);
            FibH AnnHeap(AnnHCR);
            AnnDijkstra.heap(AnnHeap, AnnHCR);
            AnnDijkstra.run(AnnStart, AnnEnd);
            if (AnnDijkstra.reached(AnnEnd))
            {
                percorso = true;
            }
        }
        mem_massimi.pop_front();
    }
    if (verbose)
    {
        save("Graph.dot");
    }
}
// fine mio codice --------------------------------------------------------------------------------------------------------------------------

std::list<std::pair<int, std::list<Mem>>> MemsGraph::visit(const SplicingGraph &sg)
{
    std::list<std::pair<int, std::list<Mem>>> paths;
    std::list<Mem> AnnPath1;
    std::list<Mem> AnnPath2;
    std::list<Mem> NovPath;
    bool FoundAnnotated = false;
    int AnnW1 = K2 + 1;
    int AnnW2 = K2 + 1;
    int NovW = K2 + 1;

    //Visiting Annotated Graph
    lemon::Dijkstra<Graph, lemon::ListDigraph::ArcMap<int>>::SetStandardHeap<FibH>::SetHeap<FibH, FibM>::Create AnnDijkstra(AnnGraph, AnnEdgesMap);
    FibM AnnHCR(AnnGraph);
    FibH AnnHeap(AnnHCR);
    AnnDijkstra.heap(AnnHeap, AnnHCR);
    AnnDijkstra.run(AnnStart, AnnEnd);
    if (AnnDijkstra.reached(AnnEnd))
    {
        AnnW1 = AnnDijkstra.dist(AnnEnd);
        if (AnnW1 <= K2)
        {
            FoundAnnotated = true;
            Path p = AnnDijkstra.path(AnnEnd);
            bool first = true;
            for (Path::ArcIt it(p); it != lemon::INVALID; ++it)
            {
                if (first)
                {
                    AnnGraph.erase(it);
                    first = false;
                }
                Arc e = it;
                Node target = AnnGraph.target(e);
                Mem m = AnnNodesMap[target];
                if (AnnGraph.id(target) != AnnGraph.id(AnnEnd))
                {
                    AnnPath1.push_back(m);
                }
            }
            paths.push_back(std::make_pair(AnnW1, AnnPath1));
        }
    }

    AnnDijkstra.run(AnnStart, AnnEnd);
    if (AnnDijkstra.reached(AnnEnd))
    {
        AnnW2 = AnnDijkstra.dist(AnnEnd);
        if (AnnW2 <= K2)
        {
            Path p = AnnDijkstra.path(AnnEnd);
            for (Path::ArcIt it(p); it != lemon::INVALID; ++it)
            {
                Arc e = it;
                Node target = AnnGraph.target(e);
                Mem m = AnnNodesMap[target];
                if (AnnGraph.id(target) != AnnGraph.id(AnnEnd))
                {
                    AnnPath2.push_back(m);
                }
            }
            paths.push_back(std::make_pair(AnnW2, AnnPath2));
        }
    }

    //Visiting Novel Graph
    lemon::Dijkstra<Graph, lemon::ListDigraph::ArcMap<int>>::SetStandardHeap<FibH>::SetHeap<FibH, FibM>::Create NovDijkstra(NovGraph, NovEdgesMap);
    FibM NovHCR(NovGraph);
    FibH NovHeap(NovHCR);
    NovDijkstra.heap(NovHeap, NovHCR);

    NovDijkstra.run(NovStart, NovEnd);
    if (NovDijkstra.reached(NovEnd))
    {
        NovW = NovDijkstra.dist(NovEnd);
        if (!FoundAnnotated)
        {
            if (NovW <= K2)
            {
                Path p = NovDijkstra.path(NovEnd);
                for (Path::ArcIt it(p); it != lemon::INVALID; ++it)
                {
                    Arc e = it;
                    Node target = NovGraph.target(e);
                    Mem m = NovNodesMap[target];
                    if (NovGraph.id(target) != NovGraph.id(NovEnd))
                    {
                        NovPath.push_back(m);
                    }
                }
            }
            if (paths.size() > 1)
                paths.pop_back();
            paths.push_front(std::make_pair(NovW, NovPath));
        }
    }
    return paths;
}

void MemsGraph::save(const std::string &s)
{
    std::ofstream myfile;

    myfile.open("Ann" + s);

    std::string dot = "digraph G {\n graph [splines=true overlap=false]\n node  [shape=ellipse, width=0.3, height=0.3]\n";
    for (NodeIt n(AnnGraph); n != lemon::INVALID; ++n)
    {
        dot += " " + std::to_string(AnnGraph.id(n)) + " [label=\"" + AnnNodesMap[n].toStr() + "\"];\n";
    }
    for (ArcIt a(AnnGraph); a != lemon::INVALID; ++a)
    {
        dot += " " + std::to_string(AnnGraph.id(AnnGraph.source(a))) + " -> " + std::to_string(AnnGraph.id(AnnGraph.target(a))) + "[label=\"" + std::to_string(AnnEdgesMap[a]) + "\"];\n";
    }
    dot += "}";

    myfile << dot;
    myfile.close();

    myfile.open("Nov" + s);

    dot = "digraph G {\n graph [splines=true overlap=false]\n node  [shape=ellipse, width=0.3, height=0.3]\n";
    for (NodeIt n(NovGraph); n != lemon::INVALID; ++n)
    {
        dot += " " + std::to_string(NovGraph.id(n)) + " [label=\"" + NovNodesMap[n].toStr() + "\"];\n";
    }
    for (ArcIt a(NovGraph); a != lemon::INVALID; ++a)
    {
        dot += " " + std::to_string(NovGraph.id(NovGraph.source(a))) + " -> " + std::to_string(NovGraph.id(NovGraph.target(a))) + "[label=\"" + std::to_string(NovEdgesMap[a]) + "\"];\n";
    }
    dot += "}";

    myfile << dot;
    myfile.close();
}