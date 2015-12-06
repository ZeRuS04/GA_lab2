#include "GeneticAlgorithm.h"
#include <qmath.h>
#include <QDebug>
#include <QFile>
#include <QLineF>
#include <QRegExp>

double fRand(double fMin, double fMax)
{
    double f = (double)rand() / RAND_MAX;
    return fMin + f * (fMax - fMin);
}

double func(double x, double y){
    return (-x*qSin(qSqrt(fabs(x)))) + (-y*qSin(qSqrt(fabs(y))));
}

GeneticAlgorithm::GeneticAlgorithm()
    : m_crossingType(MINMAX_CROSSOVER)
    , m_selectionType(TOURNEY)
    , m_count(0)
{
    m_data.resize(GA_POWER * 3 * 4);
}

void GeneticAlgorithm::fitnesFunction(gene* gene)
{
    gene->fitness = 0;
    double sum = 0;
    for(int i = 0; i < gene->length; i++)
        sum += (-gene->alleles[i]*qSin(qSqrt(fabs(gene->alleles[i]))));
    gene->fitness = sum;
}

void GeneticAlgorithm::crossOver(gene* parent1, gene* parent2)
{
    QVector<gene> g;
    g.append(gene());
    g.append(gene());
    g.append(gene());
    g.append(gene());
    for(int i = 0; i < 4; i++) {
        g[i].alleles = new double[GA_N];
        g[i].length = GA_N;
    }
    double p = fRand(0, 1);
    switch(m_crossingType) {
    case MINMAX_CROSSOVER:
        for(int i = 0; i < GA_N; i++) {
            g[0].alleles[i] = p * parent1->alleles[i] + (1-p) * parent2->alleles[i];
        }
        for(int i = 0; i < GA_N; i++) {
            g[1].alleles[i] = p * parent2->alleles[i] + (1-p) * parent1->alleles[i];
        }
        for(int i = 0; i < GA_N; i++) {
            g[2].alleles[i] = qMax(parent1->alleles[i], parent2->alleles[i]);
        }
        for(int i = 0; i < GA_N; i++) {
            g[3].alleles[i] = qMin(parent1->alleles[i], parent2->alleles[i]);
        }

        for(int i = 0; i < 4; i++) {
            fitnesFunction(&g[i]);
        }
        qSort(g);
        m_newGens << g[0] << g[1];
        delete g[2].alleles;
        delete g[3].alleles;
        break;
    }
}

void GeneticAlgorithm::initGenerator()
{
    m_genotype.clear();
    for(int i = 0; i <  GA_POWER; i++) {
        gene g;
        g.alleles = new double[GA_N];
        g.length = GA_N;
        QVector<double> usedNums;
        for (int i = 0; i < GA_N; i++) {
            g.alleles[i] = fRand(-500.,500.);
        }
        fitnesFunction(&g);
        m_genotype << g;
    }
}

void GeneticAlgorithm::selection()
{
    int parentsCount = 0;
    gene parentsArray[2];
    m_newGens.clear();

    switch (m_selectionType) {
    case ROULETTE_WHEEL:{
        qSort(m_genotype);
        double *wheel = new double[GA_POWER];
        wheel[0] = 1/m_genotype.at(0).fitness;    //Значение ФитнессФункции для 1-ого генома
        for (int i = 1; i < GA_POWER; i++){
            wheel[i] = wheel[i-1] + 1/m_genotype.at(i).fitness;   //Значение ФитнессФункции для i-ого генома
        }
        double all = wheel[GA_POWER-1];

        for (int i = 0; i < GA_POWER; i++){
            double chance = fRand(0,1);
            if(chance > GA_P_CROSS)
                continue;

            double index = fRand(0,1) * all;
            int l = 0;
            int r = GA_POWER-1;
            int c = 0;
            while (l < r){
                if(r - l == 1)
                    break;
                c = (l+r) >> 1;
                if (wheel[c] < index)
                    l = c;
                else
                    r = c;
            }

            parentsArray[parentsCount] = m_genotype.at(l);
            parentsCount++;
            if(parentsCount == 2) {
                crossOver(&parentsArray[0],&parentsArray[1]);
                parentsCount = 0;
            }
        }
    }break;
    case TOURNEY: {
        QVector<gene> genotype = m_genotype;
        while (genotype.length() >= 2) {
            int index1 = rand()%(genotype.length() - 1);
            gene g1 = genotype.at(index1);
            genotype.remove(index1);
            int index2;
            if(genotype.length() == 1)
                index2 = 0;
            else
                index2 = rand()%(genotype.length() - 1);
            gene g2 = genotype.at(index2);
            genotype.remove(index2);

            double fr1 = g1.fitness; //Значение ФитнессФункции для index1 Генома
            double fr2 = g2.fitness; //Значение ФитнессФункции для index2 Генома

            parentsArray[parentsCount] = fr1 < fr2 ? g1
                                                   : g2;
            parentsCount++;

            if(parentsCount == 2) {
                crossOver(&parentsArray[0],&parentsArray[1]);
                parentsCount = 0;
            }
        }
    } break;
    }
}

void GeneticAlgorithm::reductionOperator()
{
    m_genotype << m_newGens;
    qSort(m_genotype);
    while(m_genotype.length() > GA_POWER) {
        delete m_genotype.last().alleles;
        m_genotype.removeLast();
    }
}

void GeneticAlgorithm::mutationOperator()
{
    for(int i = 0; i < m_newGens.length(); i++) {
        double chance = fRand(0,1);
        if(chance <= GA_P_MUTATE){
            int ind = rand()%(GA_N - 1);
            m_newGens[i].alleles[ind] = fRand(-500.,500.);
        }
    }
}

void GeneticAlgorithm::run()
{

    qDebug() << "POWER:" << GA_POWER;
    qDebug() << "GENERATION_COUNT:" << GA_GENERATION_COUNT;
    qDebug() << "POWER:" << GA_POWER;
    qDebug() << "P_CROSS:" << GA_P_CROSS;
    qDebug() << "P_MUTATE:" << GA_P_MUTATE;
    initGenerator();
    if(GA_N == 2) {
        m_count = 0;
        foreach(gene g, m_genotype) {
            double x = g.alleles[0];
            double y = g.alleles[1];
            double z = g.fitness;

            double x1 = x-2;
            double y1 = y-2;
            double z1 = func(x1,y1)+1;

            double x2 = x-2;
            double y2 = y+2;
            double z2 = func(x2,y2)+1;

            double x3 = x+2;
            double y3 = y+2;
            double z3 = func(x3,y3)+1;

            double x4 = x+2;
            double y4 = y-2;
            double z4 = func(x4,y4)+1;
            add(QVector3D(x1/1500,y1/1500,z1/5000));
            add(QVector3D(x2/1500,y2/1500,z2/5000));
            add(QVector3D(x3/1500,y3/1500,z3/5000));
            add(QVector3D(x4/1500,y4/1500,z4/5000));
        }
        m_history.append(m_data);
        emit update();

        sleep(10);
    }
    for (int i = 0; i < GA_GENERATION_COUNT; i++) {
        selection();
        mutationOperator();
        reductionOperator();

        if(GA_N == 2) {
            m_count = 0;
            foreach(gene g, m_genotype) {
                double x = g.alleles[0];
                double y = g.alleles[1];

                double x1 = x-2;
                double y1 = y-2;
                double z1 = func(x1,y1)+1;

                double x2 = x-2;
                double y2 = y+2;
                double z2 = func(x2,y2)+1;

                double x3 = x+2;
                double y3 = y+2;
                double z3 = func(x3,y3)+1;

                double x4 = x+2;
                double y4 = y-2;
                double z4 = func(x4,y4)+1;
                add(QVector3D(x1/1500,y1/1500,z1/5000));
                add(QVector3D(x2/1500,y2/1500,z2/5000));
                add(QVector3D(x3/1500,y3/1500,z3/5000));
                add(QVector3D(x4/1500,y4/1500,z4/5000));
            }
            m_history.append(m_data);
            emit update();
            sleep(2);
        }

        QVector<gene> tmp = m_genotype;
        tmp.removeAll(tmp.first());
        if (tmp.isEmpty()) {
            qDebug() << "EARLY END";
            qDebug() << "Total iterations: " << i;
            break;
        }
    }
    qDebug() << "Result:";
    for( int i = 0; i <  GA_N; i++)
        qDebug() << QString("coord[%1] = ").arg(i) << m_genotype.first().alleles[i];

    qDebug() << "f = " << m_genotype.first().fitness;

}

void GeneticAlgorithm::add(const QVector3D &v)
{
    GLfloat *p = m_data.data() + m_count;
    *p++ = v.x();
    *p++ = v.y();
    *p++ = v.z();
    m_count += 3;
}

QVector<gene> GeneticAlgorithm::genotype() const
{
    return m_genotype;
}

void GeneticAlgorithm::setGenotype(const QVector<gene> &genotype)
{
    m_genotype = genotype;
}
