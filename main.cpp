#include <iostream>
#include <cstring>
#include <future> 
#include <chrono>
#include <memory>
#include <fstream>
#include <iomanip>
#include <thread>
#include <cmath>
#include <map>
#include <cassert>
#include <algorithm>

using namespace std::chrono_literals;

template <class A, class B> 
class PairWIO : public std::pair<A, B> {
    public:
        PairWIO() : std::pair<A, B>() {
            this->first = 0;
            this->second = 0; 
        }
        PairWIO(A firstI, B secondI) : std::pair<A, B>() {
            this->first = firstI;
            this->second = secondI; 
        }
        friend std::ostream& operator<< (std::ostream &out, const PairWIO<A,B> &point) {
            out << std::fixed << std::setprecision(2) << '(' << point.first << ", " << point.second << ')';
            return out;
        }
        friend std::istream& operator>> (std::istream &in, PairWIO<A,B> &point) {
            in >> point.first;
            in >> point.second;
            return in;
        } 
};

template <class T>
double distance(PairWIO<T,T> one, PairWIO<T,T> two) {
    return sqrt((one.first - two.first) * (one.first - two.first) + (one.second - two.second) * (one.second - two.second));
}

template <class T> 
class Figure {
    public:
        using Point = PairWIO<int,int>;
        Point points[4]; 
};

template <class T> 
class Rectangle : public Figure<T> {
    public:
        using Point = PairWIO<int,int>;
        Rectangle(Point ldI, Point ruI) {
            Figure<T>::points[3] = ldI;
            Figure<T>::points[1] = ruI;
            Figure<T>::points[2].first = Figure<T>::points[1].first;
            Figure<T>::points[2].second = Figure<T>::points[3].second;
            Figure<T>::points[0].first = Figure<T>::points[3].first;
            Figure<T>::points[0].second = Figure<T>::points[1].second;
        }
    
};

template <class T> 
class Rhombus : public Figure<T> {
    public:
        using Point = PairWIO<T,T>;
        Rhombus(Point in1, Point in2, int length) {
            Figure<T>::points[3] = in1;
            Figure<T>::points[1] = in2;
            double tgA = abs(Figure<T>::points[3].second - Figure<T>::points[1].second) / abs(Figure<T>::points[3].first - Figure<T>::points[1].first);
            Point mid((Figure<T>::points[3].first +Figure<T>::points[1].first) / 2, (Figure<T>::points[3].second +Figure<T>::points[1].second) / 2);
            double midLen = sqrt(length*length - pow(distance(mid, in1), 2));
            if (Figure<T>::points[3].first > Figure<T>::points[1].first) {
                Figure<T>::points[0].first = mid.first + midLen * cos(90 * M_PI / 180 - atan(tgA));
                Figure<T>::points[2].first = mid.first - midLen * cos(90 * M_PI / 180 - atan(tgA));
            } else {
                Figure<T>::points[0].first = mid.first - midLen * cos(90 * M_PI / 180 - atan(tgA));
                Figure<T>::points[2].first = mid.first + midLen * cos(90 * M_PI / 180 - atan(tgA));
            }
            if (Figure<T>::points[3].second > Figure<T>::points[1].second) {
                Figure<T>::points[0].second = mid.second - midLen * sin(90 * M_PI / 180 - atan(tgA));
                Figure<T>::points[2].second = mid.first + midLen * sin(90 * M_PI / 180 - atan(tgA));
            } else {
                Figure<T>::points[0].second = mid.second + midLen * sin(90 * M_PI / 180 - atan(tgA));
                Figure<T>::points[2].second = mid.second - midLen * sin(90 * M_PI / 180 - atan(tgA));
            }
        }
};

template <class T> 
class Trap : public Figure<T> {
    public:
        using Point = PairWIO<T,T>;
        Trap(Point in1, Point in2, double angle, double length) {
            Figure<T>::points[3] = in1;
            Figure<T>::points[2] = in2;
            double A = atan(abs(Figure<T>::points[3].second - Figure<T>::points[2].second) / abs(Figure<T>::points[3].first - Figure<T>::points[2].first));
            if (Figure<T>::points[3].first < Figure<T>::points[2].first) {
                if (Figure<T>::points[3].second < Figure<T>::points[2].second) {
                    Figure<T>::points[0].first = Figure<T>::points[3].first + length * cos(angle * M_PI / 180 + A);
                    Figure<T>::points[0].second = Figure<T>::points[3].second + length * sin(angle * M_PI / 180 + A);
                    Figure<T>::points[1].first = Figure<T>::points[2].first + length * cos(M_PI + A - angle * M_PI / 180);
                    Figure<T>::points[1].second = Figure<T>::points[2].second + length * sin(M_PI + A - angle * M_PI / 180);
                } else {
                    Figure<T>::points[0].first = Figure<T>::points[3].first + length * cos(angle * M_PI / 180 - A);
                    Figure<T>::points[0].second = Figure<T>::points[3].second + length * sin(angle * M_PI / 180 - A);
                    Figure<T>::points[1].first = Figure<T>::points[2].first + length * cos(M_PI - A - angle * M_PI / 180);
                    Figure<T>::points[1].second = Figure<T>::points[2].second + length * sin(M_PI - A - angle * M_PI / 180);
                }
            } else {
                if (Figure<T>::points[3].second > Figure<T>::points[2].second) {
                    Figure<T>::points[0].first = Figure<T>::points[3].first + length * cos(M_PI + angle * M_PI / 180 + A);
                    Figure<T>::points[0].second = Figure<T>::points[3].second + length * sin(M_PI + angle * M_PI / 180 + A);
                    Figure<T>::points[1].first = Figure<T>::points[2].first + length * cos(A - angle * M_PI / 180);
                    Figure<T>::points[1].second = Figure<T>::points[2].second + length * sin(A - angle * M_PI / 180);
                } else {
                    Figure<T>::points[0].first = Figure<T>::points[3].first + length * cos(M_PI - A + angle * M_PI / 180);
                    Figure<T>::points[0].second = Figure<T>::points[3].second + length * sin(M_PI - A + angle * M_PI / 180);
                    Figure<T>::points[1].first = Figure<T>::points[2].first + length * cos(- angle * M_PI / 180 - A);
                    Figure<T>::points[1].second = Figure<T>::points[2].second + length * sin(- angle * M_PI / 180 - A);
                }
            }
        }
};

template <typename T>
constexpr bool IsTuple = false;
template<typename ... types>
constexpr bool IsTuple<std::tuple<types...>>   = true;

template <class T, 
typename  std::enable_if<std::tuple_size<T>::value == 4>::type* = nullptr> 
void printCoor(T figure) {
    std::cout << "1st = " << std::get<0>(figure) << "\t2nd = " << std::get<1>(figure) << "\n4rd = " << std::get<3>(figure) << "\t3th = " << std::get<2>(figure) << '\n';
}

template <class T, 
typename  std::enable_if<!(IsTuple<T>)>::type* = nullptr> 
std::ostream& printCoor(std::ostream& out, T figure) {
    out << "1st = " << figure.points[0] << "\t2nd = " << figure.points[1] << "\n4rd = " << figure.points[3] << "\t3th = " << figure.points[2] << '\n';
    return out;
}

template <class T, 
typename  std::enable_if<!(IsTuple<T>)>::type* = nullptr> 
std::ofstream& printCoor(std::ofstream& out, T figure) {
    out << "1st = " << figure.points[0] << "\t2nd = " << figure.points[1] << "\n4rd = " << figure.points[3] << "\t3th = " << figure.points[2] << '\n';
    return out;
}

template <class T, 
typename  std::enable_if<std::tuple_size<T>::value == 4>::type* = nullptr>
auto centr(T figure) {
    PairWIO<double,double> out;
   
    out.first += std::get<0>(figure).first;
    out.second += std::get<0>(figure).second;
    out.first += std::get<1>(figure).first;
    out.second += std::get<1>(figure).second;
    out.first += std::get<2>(figure).first;
    out.second += std::get<2>(figure).second;
    out.first += std::get<3>(figure).first;
    out.second += std::get<3>(figure).second;
    
    out.first /= 4;
    out.second /= 4;
    return out;
}

template <class T, 
typename  std::enable_if<!(IsTuple<T>)>::type* = nullptr>
auto centr(T figure) {
    PairWIO<double,double> out;
    for (int i = 0; i < 4; i++) {
        out.first += figure.points[i].first;
        out.second += figure.points[i].second;
    }
    out.first /= 4;
    out.second /= 4;
    return out;
}

template <class T>
double geron(PairWIO<T,T> one, PairWIO<T,T> two, PairWIO<T,T> three) {
    double a = distance(one, two);
    double b = distance(two, three);
    double c = distance(one, three);
    double p = (a + b + c) / 2;
    return sqrt(p * (p - a) * (p - b) * (p - c));
}

template <class T, 
typename  std::enable_if<!(IsTuple<T>)>::type* = nullptr>
double area(T figure) { 
    return geron(figure.points[0], figure.points[1], figure.points[2]) + geron(figure.points[0], figure.points[3], figure.points[2]);
}

template <class T, 
typename  std::enable_if<std::tuple_size<T>::value == 4>::type* = nullptr>
double area(T figure) {
    return geron(std::get<0>(figure), std::get<1>(figure), std::get<2>(figure)) + geron(std::get<0>(figure), std::get<3>(figure), std::get<2>(figure));
}

template <typename T>
class TVector { 
    public:
        using value_type = T;
        using iterator = value_type*;
        
        TVector(): 
            already_used_(0), storage_size_(0), storage_(nullptr)
        {
        }

        TVector(int size, const value_type& default_value = value_type()):
            TVector()
        {
            assert(size >= 0);

            if (size == 0) {
                return;
            }

            already_used_ = size;
            storage_size_ = size;
            storage_ = std::make_unique<value_type[]>(size);

            std::fill(storage_.get(), storage_.get() + already_used_, default_value);
        }

        int size() const
        {
            return already_used_;
        }

        bool empty() const
        {
            return size() == 0;
        }

        iterator begin() const
        {
            return storage_.get();
        }
        
        iterator end() const
        {
            if (storage_.get()) {
                return storage_.get() + already_used_;
            }

            return nullptr;
        }
        void insert(iterator pos, value_type val) {
            if (already_used_ < storage_size_) {
                std::copy(pos, storage_.get() + already_used_, pos + 1);
                *pos = val;
                ++already_used_;
                return;
            }
            int next_size = 1;
            if (storage_size_) {
                next_size = storage_size_ * 2;
            }
            TVector next(next_size);
            next.already_used_ = already_used_;

            if (storage_.get()) {
                std::copy(storage_.get(), storage_.get() + storage_size_, next.storage_.get());
            }
            next.insert(pos, val);
            Swap(*this, next);
        }
        
        void erase(iterator pos) {
            std::copy(pos + 1, storage_.get() + already_used_, pos);
            --already_used_;
        }
        
        friend void Swap(TVector& lhs, TVector& rhs)
        {
            using std::swap;

            swap(lhs.already_used_, rhs.already_used_);
            swap(lhs.storage_size_, rhs.storage_size_);
            swap(lhs.storage_, rhs.storage_);
        }

        TVector& operator=(TVector other)
        {
            Swap(*this, other);
            return *this;
        }

        TVector(const TVector& other):
            TVector()
        {
            TVector next(other.storage_size_);
            next.already_used_ = other.already_used_;

            if (*(other.storage_) ) {
                std::copy(other.storage_.get(), other.storage_.get() + other.storage_size_,
                        next.storage_.get());
            }

            swap(*this, next);
        }

        ~TVector()
        {
            storage_size_ = 0;
            already_used_ = 0;
        }

        void push_back(const value_type& value)
        {
            if (already_used_ < storage_size_) {
                storage_[already_used_] = value;
                ++already_used_;
                return;
            }
            int next_size = 1;
            if (storage_size_) {
                next_size = storage_size_ * 2;
            }
            TVector next(next_size);
            next.already_used_ = already_used_;

            if (storage_.get()) {
                std::copy(storage_.get(), storage_.get() + storage_size_, next.storage_.get());
            }
            next.push_back(value);
            Swap(*this, next);
        }
        
        value_type& At(int index)
        {
            if (index < 0 || index > already_used_) {
                std::cout << "\nxxxxxx\n";
                throw std::out_of_range("You are doing this wrong!");
            }

            return storage_[index];
        }

        value_type& operator[](int index)
        {
            return At(index);
        }

    private:
        int already_used_;
        int storage_size_;
        std::unique_ptr<value_type[]> storage_;
};

template <class T>
class Subscriber {
    char* out_file;
    
    public:
        Subscriber() = default;
        Subscriber(char* out_file) : out_file(out_file) {
            if (strcmp(out_file, "stdout") != 0) {
                std::ofstream fl (out_file, std::ofstream::out);
                fl.close();
            }
        }

        int Process(TVector<Figure<T>*>& figs, int size) {
            
            
            if (strcmp(out_file, "stdout") != 0) {
                std::ofstream fl (out_file, std::ofstream::app);
                fl << size << '\n';
                for (int i = 0; i < size; i++) {
                    printCoor(fl, *figs[i]);
                    fl << "Central point: " << centr(*figs[i]) << '\n';
                    fl << "Area: " << area(*figs[i]) << "\n";
                }
                
                fl.close();

            } else {
                std::cout << size << '\n';
                for (int i = 0; i < size; i++) {
                    printCoor(std::cout, *figs[i]);
                    std::cout << "Central point: " << centr(*figs[i]) << '\n';
                    std::cout << "Area: " << area(*figs[i]) << "\n";
                }
            }
            
            return true;
        }
};

template <class T>
class Publisher {
    int size;
    int max_size;
    TVector<Figure<T>*>& figs;
    TVector<Subscriber<T>> subs;

    public:
        Publisher(TVector<Figure<T>*>& figs, int max_size):
            size(0), max_size(max_size), figs(figs), subs() {}

        
        int Add(Figure<T>* fig) {
            figs.push_back(fig);   
            size++;
            if (size == max_size) {             
                for (Subscriber<T>& sub : subs) {  
                    sub.Process(figs, max_size);
                }
                this->size = 0;  

                std::this_thread::sleep_for(1000ms);
                for (int i = figs.size() - 1; i >= 0; i--) {
                    delete figs[i];
                    figs.erase(figs.begin() + i);
                }
                return true;
            }
            return false;
        }

        void Subscribe(Subscriber<T>& sub) {
            subs.push_back(sub);
        } 
};



enum { 
    ERR, ADD,
    PRINT, DEL,
    REC, EXIT,
    CENTR, AREA,
    LES_AREA,
    TRAP, RHOMB,
    SIZE, HELP
};

void wait(std::future<int>& fut) {
    std::chrono::milliseconds span (100);

    while (fut.wait_for(span)==
        std::future_status::timeout)
    std::cout << '.' << std::flush;

    int x = fut.get();
    if (x) {
        std::cout << " Done\n\n";
    } 
}

template <class T>
void printCoorFE(T In) {
    printCoor(*In);
}

void help() {
    std::cout << "\nCommands: add, quit, help\n";
    std::cout << "Supported Figures: rectangle, trap, rhombus\n";
    std::cout << "For rectangle: 2 diagonal points (0 0 2 2)\n";
    std::cout << "For trap: 2 points of founding of a trap, angle and a length of a side ( 0 0 4 4 45 2)\n";
    std::cout << "For rhombus: 2 diagonal points and a length if a side ( 0 0 0 3 5)\n";
}

int main (int argc, char* argv[]) {
    if (argc < 2) {
        std::cout << "Error: no size of buffer entered" << '\n';
        exit(0);
    }
    using T = int;
    using Point = PairWIO<int,int>;

    int size = atoi(argv[1]);

    TVector<Figure<T>*> buf;
 
    Publisher<T> pub(buf, size);

  
    char sout[] = "stdout";
    char fout[] = "out.txt";
    Subscriber<T> sub2(sout);
    Subscriber<T> sub1(fout);
    pub.Subscribe(sub1);
    pub.Subscribe(sub2);
    
    Point tmpP1, tmpP2;
    std::string comId, figType;
    int id;
    int area_key;
    double length, angle, overallArea;
    int status = 1;
    
    std::map <std::string, int> command;
    command["add"] = ADD;
    command["rec"] = REC;
    command["rectangle"] = REC;
    command["quit"] = EXIT;
    command["q"] = EXIT;
    command["trap"] = TRAP;
    command["rhomb"] = RHOMB;
    command["help"] = HELP;
    command["h"] = HELP;

    Figure<T>* fig;

    help();
    while (status) {
        
        std::cout << "Enter command: ";
        std::cin >> comId;
        std::future<int> fut;
        switch (command[comId]) {
            case ADD:
                std::cin >> figType;
                switch (command[figType]) {
                    case REC:
                        if (!( std::cin >> tmpP1 >> tmpP2)) {
                            std::cout << "Invalid Params\n";
                            break;
                        }
                        fig = new Rectangle<int>(tmpP1, tmpP2);
                        std::cout << "Rectangle added\n";
                        break;
                    case RHOMB:
                        if (!( std::cin >> tmpP1 >> tmpP2 >> length)) {
                            std::cout << "Invalid Params\n";
                            break;
                        }
                        fig = new Rhombus<int>(tmpP1, tmpP2, length);
                        std::cout << "Rhombus added\n";
                        break;
                    case TRAP:
                        if (!( std::cin >> tmpP1 >> tmpP2 >> angle >> length)) {
                            std::cout << "Invalid Params\n";
                            break;
                        }
                        fig = new Trap<int>(tmpP1, tmpP2, angle, length);
                        std::cout << "Trap added\n";
                        break;
                    case ERR:
                        std::cout << "Unknown figure\n";
                        break;
                }
                fut = std::async(&Publisher<T>::Add,&pub,fig);
                wait(fut);

                break;
            case HELP:
                help();
                break;
            case ERR:
                std::cout << "Invalid command\n";
                break;
            case EXIT:
                for (int i = 0; i < buf.size(); i++) {
                    delete buf[i];
                }
                status = 0;
                break;
        }
        while(getchar() != '\n');
        std::cin.clear();
    }
    return 0;
}
