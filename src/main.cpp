#include <iostream>
#include <vector>
#include <array>
#include <set>
#include <map>
#include <fstream>
#include <algorithm>
#include <Eigen/Dense>
#include <Eigen/Core>
#include <string>
#include <memory>



//Représentation des Points dans le maillage
struct Node{
    double x, y;
};

//Représentation des Elements = Triangles dans le maillage
struct Triangle{
    int id;
    std::array<int, 3> node_ids;        //array pour les tableaux de taille connu
    std::string marker;          //pour savoir dans quel domaine (Omega) le triangle est ; -1: pas de marker, 0: \Gamma_0, 1: \Gamma_1, etc... (plutot mettre un string plutôt qu'un int)
};

//Représentation des Facets = Segments = Arêtes dans le maillage
struct Facet{
    int id;
    std::array<int, 2> node_ids;
    std::vector<int> triangle_ids;          //associé à un triangle dans le cas où il est sur le bord et 2 triangles sinon
    bool is_boundary;
};



//Encapsulation en un Maillage
class Mesh{
private:
    std::vector<Node> nodes;
    std::vector<Triangle> triangles;
    std::vector<Facet> facets;

    std::map<std::array<int, 2>, std::vector<int>> F;     //stockage des facettes avec ses triangles associés

public:

    //constructeur qui lit directement les données du fichier et les stocke dans la structure de donnée de maillage
    // ! : format msh commence à 1 et pas à 0
    Mesh(std::string file){
        std::ifstream ifile(file, std::ios::in);
        if(ifile){
            std::string line;
            std::map<int, std::string> domain;

            while(std::getline(ifile, line)){
                //doit prendre en compte les retours chariots présents sur les lignes
                if(line == "$PhysicalNames\r"){
                    int nb_domain;
                    ifile >> nb_domain;
                    for(int i = 0; i < nb_domain; i++){
                        int dim, num;
                        std::string name;
                        ifile >> dim >> num >> name;        
                        //les guillemets sont automatiquement retiré après avoir été stocké dans une variable
                        //pas besoin de le faire manuelement
                        domain[num] = remove_quotes(name);

                    }
                }
                if(line == "$Nodes\r"){
                    int nb_nodes;
                    ifile >> nb_nodes;
                    nodes.resize(nb_nodes);
                    for(int i = 0; i < nb_nodes; i++){
                        int id;
                        double x, y, z;
                        //z est à ignorer : on travail en 2D (z = 0)
                        ifile >> id >> x >> y >> z;
                        //format msh :
                        nodes[id-1] = {x, y};
                    }
                }
                if(line == "$Elements\r"){
                    int nb_el;
                    ifile >> nb_el;
                    for(int i = 0; i < nb_el; i++){
                        int id, type, tag1, num_domain, tag2, p1, p2, p3;
                        ifile >> id >> type >> tag1 >> num_domain >> tag2 >> p1 >> p2;
                        //si on est dans le cas d'un triangle :
                        if(type == 2){
                            ifile >> p3;
                            triangles.push_back({id, {p1-1, p2-1, p3-1}, domain[num_domain]});
                            // std::cout << "Ajout du triangle avec marker: " << domain[num_domain] << std::endl;
                            addFacet(p1 - 1, p2 - 1, triangles.size() -1);
                            addFacet(p2 - 1, p3 - 1, triangles.size() -1);
                            addFacet(p3 - 1, p1 - 1, triangles.size() -1);
                        }
                    }
                }
            }
            buildConnectivity();
            ifile.close();
        }
        else{
            std::cerr<<"Cannot open file"<<std::endl;
        }
    }

    int get_nb_nodes() const {return nodes.size();}
    Eigen::Vector2d getNode(int id) const {return Eigen::Vector2d(nodes[id].x, nodes[id].y);}
    int get_nb_triangles() const {return triangles.size();}
    const Triangle &get_triangle(int i) const {return triangles[i];}


    //fonction qui supprime les guillemets que met gmsh autour du nom des domaines
    std::string remove_quotes(const std::string& str) {
        std::string result = str;
        if (!result.empty() && result.front() == '"' && result.back() == '"'){
            //supprime le premier et dernier caractère (les " ")
            result.erase(0, 1);
            result.erase(result.size() - 1);
        }
        return result;
    }

    void addFacet(int p1, int p2, int triangle_id) {
        std::array<int, 2> segment = {p1, p2};
        std::sort(segment.begin(), segment.end());  // Garantir l'ordre croissant des nœuds
    
        // Vérification si le segment existe dans F
        if (F.find(segment) != F.end()) {
            F[segment].push_back(triangle_id);  // Facette existe, ajout du triangle
        } else {
            F[segment] = {triangle_id};  // Facette inexistante, création d'une nouvelle entrée
        }
    }


    //Ajout d'un point
    void add_node(double x, double y){
        nodes.push_back({x, y});
    }

    //Ajout d'un triangle
    void add_triangle(int id, int n1, int n2, int n3, std::string domain){
        triangles.push_back({id, {n1, n2, n3}, domain});
    }

    //Construction des facets à partir des triangles
    void buildConnectivity(){
        int next_facet_id = 0;
        for(int t = 0; t < triangles.size(); t++){              //for(auto t : triangles)
            std::vector<std::array<int, 2>> segments = {
                {triangles[t].node_ids[0], triangles[t].node_ids[1]},
                {triangles[t].node_ids[1], triangles[t].node_ids[2]},
                {triangles[t].node_ids[0], triangles[t].node_ids[2]}
            };          

            for(auto segment : segments){
                if(F.find(segment) != F.end()){
                    F[segment].push_back(t);      //facet existe
                }
                else{
                    F[segment] = {t};
                }
            }
        }
        //construction du vecteur de facets
        for(auto it : F){
            std::array<int, 2> segment = it.first;
            std::vector<int> triangle_ids = it.second;
            // std::array<int, 2> node_ids = {*segment.begin(), *std::next(segment.begin())};
            bool is_boundary = (triangle_ids.size() == 1);
            facets.push_back({next_facet_id++, segment, triangle_ids, is_boundary});
        }
    }    
    
    

    //fonction qui renvoie tous les éléments du maillage dans un domaine particulier, donné en paramètre
    // ! : format msh commence à 1 et pas à 0
    std::vector<Triangle> marked_element(std::string domain){
        std::vector<Triangle> L;      //elements in domain
        for(auto t : triangles){
            if(t.marker == domain){
                L.push_back(t);
            }
        }
        return L;
    }

    std::vector<Facet> marked_facets(const std::string &marker_name){
        std::vector<Facet> res;
        for(auto f : facets){
            // Une facette est marquée si elle est sur le bord ET que le triangle associé possède le bon marqueur
            if(f.is_boundary && !f.triangle_ids.empty()){
                int triangle_id = f.triangle_ids[0];  // par convention, sur le bord : 1 seul triangle associé
                
                if(triangles[triangle_id].marker == marker_name){
                    res.push_back(f);
                }
            }
        }
        return res;
    }
   
    
    // Affichage des facettes et de leur statut
    void printFacets() const {
        std::cout << "Facettes :\n";
        for (const auto& facet : facets) {
            std::cout << "Facette " << facet.id << " : {"
                      << facet.node_ids[0] << ", " << facet.node_ids[1] << "} - "
                      << (facet.is_boundary ? "Bord" : "Interne") << "\n";
        }
    }



    //d'un point, renvoie les triangles associés
    //renvoie pas l'id car pas explicite ici -> on aurait dû ajouter l'id dans la définission d'un triangle
    std::vector<Triangle> get_Elements_for_Point(int point_id){
        std::vector<Triangle> list_tri;
        for(auto t : triangles){
            for(int p = 0; p < 3; p++){
                if(t.node_ids[p] == point_id){
                    list_tri.push_back(t);
                }
            }
        }
        return list_tri;
    }

    //meme chose en id
    std::vector<int> get_Elements_Id_for_Point(int point_id){
        std::vector<int> list_tri_id;
        for(auto t : triangles){
            for(int p = 0; p < 3; p++){
                if(t.node_ids[p] == point_id){
                    list_tri_id.push_back(t.id);
                }
            }
        }
        return list_tri_id;
    }

    
    //d'un point, renvoie l'id des facets associées
    std::vector<int> get_Facet_Id_for_Point(int point_id){
        std::vector<int> list_f_id;
        for(auto f : facets){
            for(int p = 0; p < 2; p++){
                if(f.node_ids[p] == point_id){
                    list_f_id.push_back(f.id);
                }
            }
        }
        return list_f_id;
    }

    //autre version : renvoie les facets elles-même
    std::vector<Facet> get_Facet_for_Point(int point_id){
        std::vector<Facet> list_f;
        for(auto f : facets){
            for(int p = 0; p < 2; p++){
                if(f.node_ids[p] == point_id){
                    list_f.push_back(f);
                }
            }
        }
        return list_f;
    }


    //à partir de l'id d'une facet, renvoie les triangles associés (-> id)
    std::vector<int> get_Elements_Id_for_Facet(int facet_id){
        return facets[facet_id].triangle_ids;
    }



    //à partir d'un triangle, renvoie ses 3 triangles voisins (id)
    //si un voisin est manquant (sur le bord) renvoie -1
    std::vector<int> get_Neighbor_Id_for_Triangle(int triangle_id){
        std::vector<int> neighbors;  
        Triangle triangle = triangles[triangle_id];

        std::array<int, 2> facet1 = {triangle.node_ids[0], triangle.node_ids[1]};
        std::array<int, 2> facet2 = {triangle.node_ids[0], triangle.node_ids[2]};
        std::array<int, 2> facet3 = {triangle.node_ids[2], triangle.node_ids[1]};
        std::sort(facet1.begin(), facet1.end());
        std::sort(facet2.begin(), facet2.end());
        std::sort(facet3.begin(), facet3.end());

        for(auto facet : {facet1, facet2, facet3}){
            if(F.find(facet) != F.end()){
                std::vector<int> t_facet = F[facet];   //bien l'id des triangles associé à la facet, construit avant avec buildConnectivity
                for(auto t : t_facet){
                    if(t != triangle_id){
                        neighbors.push_back(t);
                    }
                }
            }
        }
        while(neighbors.size() < 3){   //si pas de triangles trouvés (bord) -> -1
            neighbors.push_back(-1);
        }
        return neighbors;
    }





    //détect les facets sur la frontière
    //si la facet est associée à un seul triangle => sur le bord
    bool is_Facet_on_Boundary(int facet_id){
        return get_Elements_Id_for_Facet(facet_id).size() == 1;
    }


    //renvoie l'id global d'une facet à partir de son id local dans un triangle donné
    //prend en entrée un triangle ainsi que l'id local (entre 0 et 2) d'une de ses facets
    
    // int get_Facet_Id_for_Element_Local_Facet(int triangle_id, int local_facet){
    //     Triangle t = triangles[triangle_id];
    //     std::array<int, 2> facet;
    //     switch(local_facet){
    //         case 0:
    //             facet = {t.node_ids[1], t.node_ids[2]};
    //             break;
    //         case 1:
    //             facet = {t.node_ids[2], t.node_ids[0]};
    //             break;            
    //         case 2:
    //             facet = {t.node_ids[0], t.node_ids[1]};
    //             break;
    //     }
    //     std::sort(facet.begin(), facet.end());
    //     if(F.find(facet) != F.end()){
    //         return facet;  //id global de la facet //non c'ets pas ca !!!!!!
    //     }
    //     else{return -1;}
    // }  //faut retourner l'id de la facet et pas la facet elle-même

    int get_Facet_Id_for_Element_Local_Facet(int triangle_id, int local_facet){
        Triangle t = triangles[triangle_id];
        std::array<int, 2> facet;
        switch(local_facet){               //local_facet -> 0, 1, 2    pour obtenir sa facet global
            case 0:
                facet = {t.node_ids[1], t.node_ids[2]};
                break;
            case 1:
                facet = {t.node_ids[2], t.node_ids[0]};
                break;            
            case 2:
                facet = {t.node_ids[0], t.node_ids[1]};
                break;
        }
        std::sort(facet.begin(), facet.end());

        for(const auto f : facets){
            std::array<int, 2> f_nodes = f.node_ids;
            std::sort(f_nodes.begin(), f_nodes.end());
            if(f_nodes == facet){
                return f.id;
            }
            else{
                return -1;
            }
        }
    }



    //dans la class Mesh:
    //effectue la transformation géométrique entre l'élément de référence et le triangle global
    class GeometricTransformation{
    public:

        //transformation affine pour convertir un point dans le triangle de référence en un point physique dans l'élément réel
        Eigen::Vector2d map(const Eigen::Vector2d &refPoint, const std::array<Eigen::Vector2d, 3> &nodeCoords) const{
            double xi = refPoint[0];
            double eta = refPoint[1];
            return xi*nodeCoords[1] + eta*nodeCoords[2] + (1. - xi - eta)*nodeCoords[0];
        }

        //Jacobienne pour le changement de repère
        Eigen::Matrix2d jacobian(const std::array<Eigen::Vector2d, 3> &nodeCoords) const{
            Eigen::Matrix2d J;
            J.col(0) = nodeCoords[1] - nodeCoords[0]; //dx/dxi
            J.col(1) = nodeCoords[2] - nodeCoords[0]; //dx/deta
            return J;
        }
        // Eigen::Matrix2d jacobian(const std::array<Eigen::Vector2d, 3> &nodeCoords) const {
        //     Eigen::Matrix2d J;
        //     J(0, 0) = nodeCoords[1][0] - nodeCoords[0][0]; // dx/dxi
        //     J(0, 1) = nodeCoords[2][0] - nodeCoords[0][0]; // dx/deta
        //     J(1, 0) = nodeCoords[1][1] - nodeCoords[0][1]; // dy/dxi
        //     J(1, 1) = nodeCoords[2][1] - nodeCoords[0][1]; // dy/deta
        //     return J;
        // }
    };

};







//partie FEM
//classe abstraite qui représente les fonctions de base -> Lagrange
class BasisFunction{
public:
    virtual Eigen::MatrixXd evaluate(const Eigen::MatrixXd &refPoints) const = 0;   //évalue les fontions de base aux points de référence
    virtual Eigen::MatrixXd gradient(const Eigen::MatrixXd &refPoints) const = 0;   //évalue le gradient des fonctions de base aux points de référence
    virtual ~BasisFunction() = default;
};

//fonction de forme de Lagrange linéaire sur un triangle de référence
class LagrangeBasis : public BasisFunction{
public:
    //évalue les points de références dans la base de Lagrange
    Eigen::MatrixXd evaluate(const Eigen::MatrixXd &refPoints) const override{
        int N = refPoints.cols();     //nombre de points
        Eigen::MatrixXd values(3, N);    //stocke l'évaluation des 3 fonction de base sur chaque points
        values.row(0) = Eigen::MatrixXd::Ones(1, N) - refPoints.row(0) - refPoints.row(1);   //1 - x - y
        values.row(1) = refPoints.row(0);   //x
        values.row(2) = refPoints.row(1);   //y
        return values;
    }

    Eigen::MatrixXd gradient(const Eigen::MatrixXd &refPoints) const override{
        int N = refPoints.cols();
        Eigen::MatrixXd grad(3*2, N);
        grad.row(0).setConstant(-1);
        grad.row(1).setConstant(-1);
        grad.row(2).setConstant(1);
        grad.row(3).setConstant(0);
        grad.row(4).setConstant(0);
        grad.row(5).setConstant(1);
        return grad;
    }
};



//encapsule FEM function space
//définie un espace de fonction dans un maillage
class FunctionSpace{
public:
    std::shared_ptr<Mesh> mesh;
    std::unique_ptr<BasisFunction> basis;
    std::vector< std::vector<int> > dofTable;

    FunctionSpace(std::shared_ptr<Mesh> m, std::unique_ptr<BasisFunction> b) : mesh(m), basis(std::move(b)){}

    void initDofs(){
        int nbNodes = mesh->get_nb_nodes();
        dofTable.resize(mesh->get_nb_triangles());

        for(int i = 0; i < mesh->get_nb_triangles(); i++){
            auto tri = mesh->get_triangle(i);
            dofTable[i] = { tri.node_ids[0], tri.node_ids[1], tri.node_ids[2] };
        } 
    }

    //représente une fonction locale définie sur un triangle à l'aide de ses coefficients 
    class Element{
    public:
        Eigen::VectorXd weights;
        std::shared_ptr<BasisFunction> basis;

        Element(const Eigen::VectorXd &coeffs, std::shared_ptr<BasisFunction> b) : weights(coeffs), basis(std::move(b)){}

        Eigen::VectorXd evaluate(const Eigen::MatrixXd &refPoints) const{
            Eigen::MatrixXd phi = basis->evaluate(refPoints);
            return phi.transpose() *weights;
        }
    };
};




double error_L2(const FunctionSpace &space, const Eigen::VectorXd &uh, std::function<double(Eigen::Vector2d)> u_exact, const Eigen::MatrixXd &quadPoints, const Eigen::VectorXd &quadweights){
    double err_L2 = 0.;
    std::shared_ptr<Mesh> mesh = space.mesh;

    Mesh::GeometricTransformation geo;
    int Q = quadPoints.cols();

    for(int k = 0; k < mesh->get_nb_triangles(); k++){
        Triangle K = mesh->get_triangle(k);

        std::array<Eigen::Vector2d, 3> nodeCoords = {mesh->getNode(K.node_ids[0]), mesh->getNode(K.node_ids[1]), mesh->getNode(K.node_ids[2])};

        //coeffs locaux de u_h
        Eigen::VectorXd local_uh(3);
        for(int i = 0; i < 3; i++){
            local_uh(i) = uh(K.node_ids[i]);
        }

        //élément local
        FunctionSpace::Element elem(local_uh, std::make_unique<LagrangeBasis>());

        //Jacobienne
        Eigen::Matrix2d J = geo.jacobian(nodeCoords);
        double detJ = std::abs(J.determinant());


        for(int q = 0; q < Q; q++){
            Eigen::Vector2d refPoint = quadPoints.col(q);
            Eigen::Vector2d reelPoint = geo.map(refPoint, nodeCoords);

            Eigen::MatrixXd refMat(2, 1);
            refMat.col(0) = refPoint;

            double u_val = u_exact(reelPoint);
            double uh_val = elem.evaluate(refMat)(0);  //reçoit un vecteur de un élément mais nous on veux comme un double -> (0)

            err_L2 += quadweights(q)* (u_val - uh_val)*(u_val - uh_val) * detJ;
        }
    }
    return std::sqrt(err_L2);
}




double error_H1(const FunctionSpace &space, const Eigen::VectorXd &uh, std::function<double(Eigen::Vector2d)> u_exact, std::function<Eigen::Vector2d(Eigen::Vector2d)> grad_u_exact, const Eigen::MatrixXd &quadPoints, const Eigen::VectorXd &quadWeights){
    double err_grad = 0.;
    double err_L2 = error_L2(space, uh, u_exact, quadPoints, quadWeights);
    double err_L2_square = err_L2 * err_L2;

    std::shared_ptr<Mesh> mesh = space.mesh;
    Mesh::GeometricTransformation geo;

    for(int k = 0; k < mesh->get_nb_triangles(); k++){
        Triangle K = mesh->get_triangle(k);
        std::array<Eigen::Vector2d, 3> nodeCoords = {mesh->getNode(K.node_ids[0]), mesh->getNode(K.node_ids[1]), mesh->getNode(K.node_ids[2])};

        Eigen::VectorXd local_uh(3);
        for(int i = 0; i < 3; i++){
            local_uh(i) = uh(K.node_ids[i]);
        }

        FunctionSpace::Element elem(local_uh, std::make_unique<LagrangeBasis>());
        
        Eigen::Matrix2d J = geo.jacobian(nodeCoords);
        double detJ = std::abs(J.determinant());
        Eigen::Matrix2d invJt = J.inverse().transpose();

        Eigen::MatrixXd grad_ref = elem.basis->gradient(quadPoints);


        for(int q = 0; q < quadPoints.cols(); q++){
            Eigen::Vector2d refPoint = quadPoints.col(q);
            Eigen::Vector2d reelPoint = geo.map(refPoint, nodeCoords);

            Eigen::Vector2d grad_u_val = grad_u_exact(reelPoint);

            Eigen::Vector2d grad_uh_val = Eigen::Vector2d::Zero();
            for(int i = 0; i < 3; i++){
                Eigen::Vector2d grad_phi_i;
                grad_phi_i(0) = grad_ref(2*i +0, q);
                grad_phi_i(1) = grad_ref(2*i +1, q);
                grad_uh_val += local_uh(i) * invJt * grad_phi_i;
            }

            err_grad += quadWeights(q)* ((grad_u_val - grad_uh_val).dot(grad_u_val - grad_uh_val) ) *detJ;
        }
    }
    return std::sqrt(err_grad + err_L2_square);
}










// //pour chaque triangle on récupère les 3 sommets et on fait la transformation géométrique
// //calcul du jacobien : Aire K = 1/2 *detJ
// //gradient des foctions de base : nabla phi_i = J^-t * nabla * ^phi_i
// //intégration : en P1 les gradients sont constants
// void assemble(const FunctionSpace &fs, Eigen::SparseMatrix<double> &A, Eigen::VectorXd &b, std::function<double(const Eigen::Vector2d &)> f){
//     int nb_nodes = fs.mesh->get_nb_nodes;
//     A.resize(nb_nodes, nb_nodes);
//     b = Eigen::VectorXd::Zero(nb_nodes);

//     std::vector<Eigen::Triplet<double>> triplets;
//     Mesh::GeometricTransformation geo;

//     for(int e_id = 0; e_id < fs.mesh->get_nb_triangles(); e_id++){
//         auto tri = fs.mesh->get_triangle(e_id);

//         //coordonnées physiques des 3 sommets
//         std::array<Eigen::Vector2d, 3> nodeCoords = {
//             fs.mesh->get_node(tri.node_ids[0]),
//             fs.mesh->get_node(tri.node_ids[1]),
//             fs.mesh->get_node(tri.node_ids[2])
//         };

//         Eigen::Matrix2d J = geo.jacobian(nodeCoords);
//         double area = 0.5*std::abs(J.determinant());
//         Eigen::Matrix2d J_invt = J.inverse().transpose();

//         //Gradient des fonctions de base en référence :
//         Eigen::Matrix<double, 3, 2> gradref;
//         gradref << -1, -1, 1, 0, 0, 1;

//         //gradient physique:
//         Eigen::Matrix<double, 3, 2> gradphys = gradref * J_invt;

//         //intégrale locale A_ij :
//         for(int i = 0; i < 3; i++){
//             for(int j = 0; j < 3; j++){
//                 double val = area *gradphys.row(1).dot(gradphys.row(j));
//                 triplets.emplace_back(tri.node_ids[i], tri.node_ids[j], val);
//             }
//         }


//         //RHS local (source constante évaluée au centre du triangle) :
//         Eigen::Vector2d x_center = (nodeCoords[0] + nodeCoords[1] + nodeCoords[2])/3.;
//         double fval = f(x_center);
//         for(int i = 0; i < 3; i++){
//             b[tri.nodes_ids[i]] += area * fval /3. ;
//         }
//     }
//     A.setFromTriplets(triplets.begin(), triplets.end());
// }





// void compute_L2_H1_error(const Mesh &mesh, const Eigen::VectorXd &uh, std::function<double(const Eigen::Vector2d&)> u_exact,
//                     std::function<Eigen::Vector2d(const Eigen::Vector2d&)> grad_u_exact, double L2_error, double H1_error){
//     L2_error = 0.0;
//     HA_error = 0.0;

//     Mesh::GeometricTransformation geo;;
//     LagrangeBasis basis;

//     //points de quadrature au centre du triangle pour P1:
//     Eigen::MatrixXd refpoint(2, 1);
//     refpoint << 1./3., 1./3.;

//     Eigen::MatrixXd phi = basis.evaluate(refpoint);
//     Eigen::Matrix<double, 3, 2> gradref;
//     gradref << -1, -1, 1, 0, 0, 1;

//     for(int t = 0; t < mesh.get_nb_triangles(); t++){
//         const Triangle &tri = mesh.get_triangle(t);
//         std::array<Eigen::Vector2d, 3> nodeCoords {
//             mesh.get_node(tri.node_ids[0]),
//             mesh.get_node(tri.node_ids[1]),
//             mesh.get_node(tri.nodes_ids[2])
//         };

//         //uh évalué au centre :
//         double uh_val = 0.;
//         for(int i = 0; i < 3; i++){
//             uh_vak += uh[tri.node_ids[i]] * phi(i, 0);
//         }

//         //u exact:
//         Eigen::Vector2d x = (nodeCoords[0] + nodeCoords[1] + nodeCoords[2])/3.;
//         double u_val = u_exact(x);
//         Eigen::Vector2d grad_u = grad_u_exact(x);

//         //Aire :
//         Eigen::Matrix2d J = geo.jacobian(nodeCoords);
//         double area = 0.5 *std::abs(J.determinant());
//         Eigen::Matrix2d J_invt = J.inverse().transpose();
//         Eigen::Matrix<double, 3, 2> gradphys = gradref *J_invt;

//         //gradient de uh :
//         Eigen::Vector2d grad_uh = Eigen::Vector2d::Zeor();
//         for(int i = 0; i < 3; i++){
//             grad_uh += uh[tri.nodes_ids[i]] * gradphys.row(i).transpose();
//         }

//         L2_error += area*std::pow(u_val - uhval, 2);
//         H1_error += area*(grad_u - grad_uh).squaredNorm();
//     }

//     L2_error = std::sqrt(L2_error);
//     H1_error = std::sqrt(H1_error);
// }





int main() {
    // Mesh mesh;

    // mesh.add_node(0.0, 0.0);
    // mesh.add_node(1.0, 0.0);
    // mesh.add_node(0.0, 1.0);
    // mesh.add_node(1.0, 1.0);

    // mesh.add_triangle(0, 1, 2, "Omega1");
    // mesh.add_triangle(1, 3, 2, "Omega2");

    // mesh.buildConnectivity();

    // mesh.printFacets();
    

    // auto elements_omega1 = mesh.marked_element("Omega1");

    // std::cout << "Triangles dans Omega1:\n";
    // for (auto t : elements_omega1) {
    //     std::cout << "Triangle avec noeuds {"
    //                 << t.node_ids[0] << ", " << t.node_ids[1] << ", " << t.node_ids[2] 
    //                 << "} appartient à " << t.marker << "\n";
    // }
    // std::cout<<std::endl;

    Mesh mesh("../meshes/meshes_with_domain/mesh_test0.msh");
    mesh.printFacets();
    auto el = mesh.marked_element("Omega");
    std::cout << "Triangles dans Omega:\n";
    for (auto t : el) {
        std::cout << "Triangle avec noeuds {"
                    << t.node_ids[0] << ", " << t.node_ids[1] << ", " << t.node_ids[2] 
                    << "} appartient à " << t.marker << "\n";
    }
    std::cout<<std::endl;

    return 0;
}
