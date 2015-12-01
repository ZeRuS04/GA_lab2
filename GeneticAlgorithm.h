#ifndef GENETICALGORITHM_H
#define GENETICALGORITHM_H

#include <qopengl.h>
#include <QVector3D>
#include <QThread>
#include <QVector>
#include <QList>
#include <QPointF>
#include <float.h>
#include <qmath.h>

#define GA_POWER        50
#define GA_P_CROSS      0.55
#define GA_P_MUTATE     0.001
#define GA_GENERATION_COUNT     100

#define GA_N 2
enum SelectionType {
    TOURNEY,
    ROULETTE_WHEEL
};

enum CrossingType {
    MINMAX_CROSSOVER
};

struct gene {
    double *alleles;
    unsigned short length;
    double fitness;

    bool operator<(gene const & b) const {return ((this->fitness) < (b.fitness));}
    bool operator==(gene const & gn) const  {
        if(this->fitness != gn.fitness)
            return false;

        if(this->length != gn.length)
            return false;

        for(int i = 0; i < this->length; i++) {
            if(this->alleles[i] != gn.alleles[i])
                return false;
        }
        return true;
    }
};

class  GeneticAlgorithm : public QThread
{
    Q_OBJECT

    void fitnesFunction(gene* gene);
    void crossOver(gene* parent1, gene* parent2);
public:
    explicit GeneticAlgorithm();

    void initGenerator();
    void selection();
    void reductionOperator();
    void mutationOperator();

    QVector<gene> genotype() const;
    void setGenotype(const QVector<gene> &genotype);


    const GLfloat *constData() const { return m_data.constData(); }
    int count() const { return m_count; }
    int vertexCount() const { return m_count / 3; }


protected:
    // QThread interface
    void run();

signals:
    void update();

private:
    QVector<gene> m_genotype;
    QVector<gene> m_newGens;

    QVector<GLfloat> m_data;
    int m_count;

    CrossingType m_crossingType;
    SelectionType m_selectionType;

    void add(const QVector3D &v);
};

#endif // GENETICALGORITHM_H
