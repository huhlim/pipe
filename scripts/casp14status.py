#!/usr/bin/env python

import os
import sys
import cgi
import pymysql
import datetime

import base64

PRED_s = ['FEIG-S', 'FEIG-R1', 'FEIG-R2', 'FEIG-R3', 'FEIG']
STATUS_s = ['', '_R1', '_R2', '_R3', '_H']

def print_header():
    print ('''Content-type: text/html

<html>
<head>
    <title>FeigLab CASP14 status</title>
    <meta http-equiv="refresh" content="600">
</head>
''')

def get_targets_from_db():
    db = pymysql.connect(host='blue7', user='casp14', password='CA9ee04cf01SP14!', db='casp14',
                cursorclass=pymysql.cursors.DictCursor)
    with db.cursor() as cursor:
	sql = 'SELECT * from targets'
	cursor.execute(sql)
        result_s = cursor.fetchall()
    db.close()
    return result_s

def print_targets(target_s):
    now = datetime.datetime.now()
    #
    print ("<table border=1>")
    print ("  <thead>")
    print ("    <tr>")
    print ("      <th rowspan=2 width=120>Target</th>")
    print ("      <th rowspan=2 width=120>Type</th>")
    print ("      <th rowspan=2 width=120>Entry date</th>")
    print ("      <th rowspan=2 width=120>Server exp.</th>")
    print ("      <th rowspan=2 width=120>Human exp.</th>")
    print ("      <th colspan=5 width=500>Status / Time left</th>")
    print ("    </tr>")
    print ("    <tr>")
    for pred in PRED_s:
        print ("      <th width=100>%s</th>"%(pred))
    print ("    </tr>")
    print ("  </thead>")
    print ('  <tbody style="border:1px solid black;border-collaspse:collapse;border:1px;">')

    for target in sorted(target_s, reverse=True, key=lambda x: x['submit_time']):
        target_type = target['target_type']
        if target_type is None:
            target_type = 'NONE'
        submit_time = target['submit_time']
        #
        server_exp = submit_time + datetime.timedelta(days=3)
        server_exp = server_exp.replace(hour=15, minute=0, second=0, microsecond=0)
        server_left = (server_exp - now).total_seconds() / 3600.
        #
        human_exp = submit_time + datetime.timedelta(days=21)
        human_exp = human_exp.replace(hour=15, minute=0, second=0, microsecond=0)
        human_left = (human_exp - now).total_seconds() / 3600.
        #
        print ("    <tr>")
        if target['target_id'].startswith("R"):
            print ("      <td align='center'><b>%s</b></td>"%(target['target_id']))
        else:
            print ("      <td align='center'>%s</td>"%(target['target_id']))
        print ("      <td align='center'>%s</td>"%(target_type.lower()))
        print ("      <td align='center'>%s</td>"%(submit_time.strftime("%Y-%m-%d")))
        if target_type in ['HUMAN', 'SERVER', 'REFINE']:
            print ("      <td align='center'>%s</td>"%(server_exp.strftime("%Y-%m-%d")))
        else:
            print ("      <td align='center'>-</td>")
        if target_type in ['HUMAN', 'REFINE']:
            print ("      <td align='center'>%s</td>"%(human_exp.strftime("%Y-%m-%d")))
        else:
            print ("      <td align='center'>-</td>")
        for pred,status_name in zip(PRED_s,STATUS_s):
            status = target['status%s'%status_name]
            if status == 'DONE':
                if pred == 'FEIG-S' and (server_left < 0 or target_type not in ['HUMAN', 'SERVER']):
                    png_fn = '/green/s2/huhlim/work/casp14/%s/%s/trRosetta/%s.trRosetta.png'%(pred, target['target_id'], target['target_id'])
                    if os.path.exists(png_fn):
                        print ("      <td align='center' bgcolor='e0e0e0'><a href='./casp14status.py?model=%s'>O</a></td>"%target['target_id'])
                    else:
                        print ("      <td align='center' bgcolor='e0e0e0'>O</td>")
                else:
                    print ("      <td align='center' bgcolor='e0e0e0'>O</td>")
            elif (target['target_id'].startswith("T") and pred != 'FEIG-S' and target_type != 'HUMAN') or \
                    (target['target_id'].startswith("R") and pred in ['FEIG-R1', 'FEIG-R2', 'FEIG-R3', 'FEIG-R4']):
                print ("      <td align='center' bgcolor='e0e0e0'>.</td>")
            else:
                if pred == 'FEIG-S':
                    time_left = server_left
                else:
                    time_left = human_left
                status_dat = '/green/s2/huhlim/work/casp14/%s/%s/status.dat'%(pred, target['target_id'])
                #
                if time_left < 72.:
                    out = '%d h'%time_left
                else:
                    out = '%d d'%(time_left/24.)
                if status == 'RUN':
                    if os.path.exists(status_dat):
                        out = '<a href="./casp14status.py?id=%s&pred=%s"><b>%s</b></a>'%(target['target_id'], pred, out)
                    else:
                        out = '<b>%s</b>'%out
                if time_left < 0:
                    color = 'ff0000'
                elif time_left < 24.:
                    color = 'ffcc99'
                elif time_left < 72.:
                    color = 'ccff99'
                else:
                    color = '99ffff'

                print ("      <td align='center' bgcolor='%s'>%s</td>"%(color, out))
        print ("    <tr>")
    print ("  </tbody>")
    print ("</table>")

def print_status(form):
    id = form.getvalue("id")
    pred = form.getvalue("pred")
    #
    status_dat = '/green/s2/huhlim/work/casp14/%s/%s/status.dat'%(pred, id)
    with open(status_dat) as fp:
        for line in fp:
            line = line.rstrip().replace(" ",'&nbsp')
            if 'RUN' in line:
                line = "<span style='font-size:0.9em; font-weight: bold; font-family: Monospace ; color: blue ;'>%s</span>"%line
            elif 'WAIT' in line:
                line = "<span style='font-size:0.9em; font-weight: bold; font-family: Monospace ; color: black ;'>%s</span>"%line
            else:
                line = "<span style='font-family: Monospace;'>%s</span>"%line
            print ("%s</br>"%line)

def print_model(form, target_s):
    id = form.getvalue("model")
    #
    now = datetime.datetime.now()
    #
    found = False
    for target in sorted(target_s, reverse=True, key=lambda x: x['submit_time']):
        if target['target_id'] == id:
            found = True ; break
    #
    if found:
        submit_time = target['submit_time']
        #
        server_exp = submit_time + datetime.timedelta(days=3)
        server_exp = server_exp.replace(hour=15, minute=0, second=0, microsecond=0)
        server_left = (server_exp - now).total_seconds() / 3600.
        #
        png_fn = '/green/s2/huhlim/work/casp14/%s/%s/trRosetta/%s.trRosetta.png'%('FEIG-S', id, id)
        if (server_left < 0 or target['target_type'] not in ['HUMAN', 'SERVER']) and os.path.exists(png_fn):
            with open(png_fn, 'rb') as fp:
                image = fp.read()
            data = base64.b64encode(image).decode("ascii")
            #
            print ("<img src='data:image/png;base64,%s'/>"%data)
        else:
            print ("<b>Not ready to show</b>")
    else:
        print ("<b>Cannot find %s</b>"%id)

def main():
    print_header()
    #
    form = cgi.FieldStorage()
    #
    print ("<body>")
    if form.has_key("id") and form.has_key("pred"):
        print_status(form)
    elif form.has_key("model"):
        target_s = get_targets_from_db()
        print_model(form, target_s)
    else:
        target_s = get_targets_from_db()
        print_targets(target_s)
    #
    print ("</body>")
    print ("</html>")


if __name__ == '__main__':
    main()
