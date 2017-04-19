import numpy as np
import pandas as pd

def birthplace_fixer(code):
    '''
    Maps code to birthplace facility type.
    '''
    dicto = {1:'Hospital',2:'Freestanding Birth Center',3:'Home',4:'Home',
             5:'Home',6:'Clinic/Doctor\'s Office',7:'Other',9:'Unknown'}

    if not code or int(code) not in dicto.keys():
        return 'Unknown'
    else:
        return dicto[int(code)]


def race_fixer(code):
    '''
    Maps code to race.
    '''
    dicto = {1:'White',2:'Black',3:'Native American',4:'Asian',
             5:'Native Hawaiian or Pacific Islander',6:'2+ races',
             9:'Unknown'}

    if int(code) not in dicto.keys():
        return 'Unknown'
    else:
        return dicto[int(code)]


def edu_fixer(code):
    '''
    Maps code to education level.
    '''

    dicto = {1:'8th grade or less',2:'High school, no diploma',3:'HSD/GED',
             4:'Some college, no degree',5:'Associate\'s',6:'Bachelor\'s',
             7:'Master\'s',8:'Doctorate/Professional',9:'Unknown'}

    if not code:
        return 'Unknown'
    else:
        return dicto[int(code)]


def priors(integer):
    '''
    Recodes 99 as Unknown.
    '''
    if int(integer) == 99:
        return np.NaN
    else:
        return int(integer)


def str_int(string):
    '''
    Recodes 888 to First and 999 to Unknown.
    '''

    if string == '888':
        return np.NaN   #'First'
    elif string == '999' or not string:
        return np.NaN   #'Unknown'
    else:
        return int(string)


def attendant(code):
    '''
    Maps code to attendant.
    '''
    dicto = {1:'MD',2:'DO',3:'CNM',4:'Other midwife',5:'Other',9:'Unknown'}

    return dicto[int(code)]


def pay(code):
    '''
    Maps code to payment type.
    '''
    dicto = {1:'Medicaid',2:'Private Insurance',3:'Self-pay',4:'Other',
             9:'Unknown'}

    if not code:
        return 'Unknown'
    else:
        return dicto[int(code)]


def pnc_fixer1(string):
    '''
    Recodes 99 and 0 to np.NaN.
    '''

    if string == '99' or string == '0' or not string:
        return np.NaN   #'Unknown'
    else:
        return int(string)

def pnc_fixer2(string):
    '''
    Recodes 99 to np.NaN.
    '''

    if string == '99' or not string:
        return np.NaN   #'Unknown'
    else:
        return int(string)


def months_between(string):
    '''
    Drops plural deliveries.
    '''

    drop = ['000','001','002','003','888','999','3']

    if string in drop or not string:
        return np.NaN
    else:
        return int(string)


if __name__ == '__main__':

    fname = 'natl2015.csv'

    uc = ['dob_mm', 'bfacil', 'mager', 'mrace6', 'meduc', 'fagecomb', 'frace6',
          'feduc', 'priorlive', 'priordead', 'priorterm', 'illb_r', 'precare',
          'previs', 'attend', 'pay_rec', 'sex', 'dplural']

    cv = {'bfacil':birthplace_fixer,'mrace6':race_fixer,'meduc':edu_fixer,
          'frace6':race_fixer,'priorlive':priors,'priordead':priors,
          'priorterm':priors,'illb_r':months_between,'previs':str_int,
          'fagecomb':str_int,'attend':attendant,'pay_rec':pay,'feduc':edu_fixer,
          'previs':pnc_fixer2,'precare':pnc_fixer1}

    df = pd.read_csv(fname, usecols = uc, converters = cv)

    cnames = {'dob_mm':'BirthMonth','bfacil':'FacilityType',
              'mager':'MaternalAge','mrace6':'MaternalRace',
              'meduc':'MaternalEducation','fagecomb':'PaternalAge',
              'frace6':'PaternalAge','feduc':'PaternalEducation',
              'priorlive':'PriorLiveBirths','priordead':'PriorNonLiveBirths',
              'priorterm':'PriorTerminations','illb_r':'MonthsBtwnPregnancies',
              'precare':'PNCMonth','previs':'PNCVisits','attend':'Attendant',
              'pay_rec':'PaymentMethod','sex':'ChildSex','dplural':'Multiples'}

    df = df.rename(index=str, columns = cnames)

    df = df.drop(['ChildSex','PaymentMethod','Attendant',
                  'PaternalEducation','PaternalAge','PaternalRace',
                  'MaternalRace','BirthMonth','FacilityType'], axis=1)

    df_a_l = df[df.MaternalEducation != 'Unknown']
    df_a_l = df_a_l[df_a_l.MaternalEducation != 'Doctorate/Professional']
    df_a_l = df_a_l[df_a_l.MaternalEducation != 'Master\'s']
    df_a_l = df_a_l[df_a_l.MaternalEducation != 'Bachelor\'s']

    df_b_h = df[df.MaternalEducation != 'Unknown']
    df_b_h = df_b_h[df_b_h.MaternalEducation != '8th grade or less']
    df_b_h = df_b_h[df_b_h.MaternalEducation != 'High school, no diploma']
    df_b_h = df_b_h[df_b_h.MaternalEducation != 'HSD/GED']
    df_b_h = df_b_h[df_b_h.MaternalEducation != 'Some college, no degree']
    df_b_h = df_b_h[df_b_h.MaternalEducation != 'Associate\'s']


    df.to_csv('df.csv', index=False)
    df_a_l.to_csv('dfal.csv', index=False)
    df_b_h.to_csv('dfbh.csv', index=False)
