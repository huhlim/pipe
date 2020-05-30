#!/usr/bin/env python

import os
import sys
import cgi
import pymysql

MAX_RUNNING_JOB = 10 

def print_header():
    print ('''Content-type: text/html

<html>
    <head>
        <title>FeigLab CASP14 submission</title>
    </head>
''')

def print_body_error():
    print ("<body>")
    print ("    <a><b>At least one of the required field is missing.</b></a>")
    print ("</body>")

def print_body(form, status):
    print ("<body>")
    if status == 0:
        print ("    <a><b>A CASP14 target has been submitted to FeigLab server.</b></a></br>")
        print ("    <a>Target ID: %s</a></br>"%form.getvalue("TARGET"))
        print ("    <a>EMAIL: %s</a></br>"%form.getvalue("EMAIL"))
        print ("    <a>SEQUENCE: %s</a></br>"%form.sequence)
    elif status == 1:
        print ("    <a><b>Invalid target.</b></a></br>")
    elif status == 2:
        print ("    <a><b>The submitted target and the corresponding sequence is already in the database.</b></a></br>")
    elif status == 3:
        print ("    <a><b>You have submitted too many jobs.</b></a></br>")
    print ("</body>")

def submit_to_db(form):
    target_id = form.getvalue("TARGET")
    email = form.getvalue("EMAIL")
    if target_id[:2] not in ['R0','R1','T0','T1','T9'] and email != 'huhlim@gmail.com':
        return 1
    #
    # assume that there will be only one protein in the SEQUENCE
    sequence = []
    for line in form.getvalue("SEQUENCE").split("\n"):
        if line.startswith(">"): continue
        sequence.append(line.strip())
    sequence = ''.join(sequence)
    form.sequence = sequence
    ip_addr = cgi.escape(os.environ['REMOTE_ADDR'])
    #
    db = pymysql.connect(host='blue7', user='casp14', password='CA9ee04cf01SP14!', db='casp14',
                cursorclass=pymysql.cursors.DictCursor)
    #
    with db.cursor() as cursor:
        sql = 'SELECT * from targets'
        cursor.execute(sql)
        result_s = cursor.fetchall()
        #
        is_new = True
        job_id = 0
        n_running = 0
        for result in result_s:
            if result['target_id'] == target_id and result['sequence'] == sequence and result['email'] == email:
                is_new = False
                break
            else:
                job_id = max(job_id, result['id']+1)
                if result['ip_addr'] == ip_addr and result['status'] in ['SUBMIT', 'RUN']:
                    n_running += 1
    #
    if not is_new:
        return 2
    elif n_running+1 > MAX_RUNNING_JOB:
        return 3
    #
    with db.cursor() as cursor:
        sql = 'INSERT INTO targets (id, target_id, sequence, submit_time, email, ip_addr, status) '
        sql += 'VALUES (%d, "%s", "%s", NOW(), "%s", "%s", "SUBMIT")'%(job_id, target_id, sequence, email, ip_addr)
        cursor.execute(sql)
    db.commit()
    #
    db.close()
    return 0

def print_submit_form():
    print ("""
    <body>
        <form name='submit_form' action='./casp14submit' method='post'>
            <table>
                <tr>
                    <td>Target Name</td>
                    <td>
                        <input type='text' name='TARGET' value=""/>
                    </td>
                </tr>
                <tr>
                    <td>E-mail</td>
                    <td>
                        <input type='text' name='EMAIL' value=""/>
                    </td>
                </tr>
                <tr>
                    <td>Sequence</td>
                    <td>
                        <input type='text' name='SEQUENCE' value=""/>
                    </td>
                </tr>
                <tr>
                    <td>
                        <input type='submit' name='submit' value='submit' />
                        <input type="reset" name="reset" value="reset" />
                    </td>
                </tr>
            </table>
        </form>
    </body>""")

def main():
    print_header()
    #
    form = cgi.FieldStorage()
    #print (form.keys())
    if not form.has_key("TARGET"):
        submit_status = False
    elif not form.has_key("SEQUENCE"):
        submit_status = False
    elif not form.has_key("EMAIL"):
        submit_status = False
    else:
        submit_status = True
    #
    if submit_status:
        status = submit_to_db(form)
        print_body(form, status)
    else:
        #print_body_error()
        print_submit_form()
    #
    print ("</html>")


if __name__ == '__main__':
    main()
