""" Server for Maxwell.

    Receives user jobs (simulations), executes them, and returns the results.
"""

import BaseHTTPServer
import cgi, sys, os, shutil, uuid, time
import subprocess, shlex
from SocketServer import ForkingMixIn

PAUSE_TIME = 0.1 # Polling interval for the simulation status file.

project_template = """
name %s
oticket 0
fshare %d
acl NONE
xacl NONE
"""

PASSWORD = sys.argv[2]

class MaxwellHandler(BaseHTTPServer.BaseHTTPRequestHandler):
    """ Handler for the server.

        Interacts with the client in the following way:
        *   Client POSTs simulation input data.
        *   Simulation is submitted to the scheduler (sge).
        *   Simulation is executed.
        *   Client is redirected to new url for GET.
        *   Client GETs simulation results.
    """

    def do_POST(self):
        """ Handles input simulation data.

            Basic operation involves:
            *   Check client credentials.
            *   Save input data to file.
            *   Submit job to scheduler.
            *   Stream status reports to client.
            *   On job completion, return output files to client.

            Should handle the following cases:
            *   If client disconnects, kill job.
            *   Job duration time-limit? 
        """
        # Interpret the structure of the post message.
        form = cgi.FieldStorage( \
            fp=self.rfile, \
            headers=self.headers, \
            environ={'REQUEST_METHOD':'POST', \
                    'CONTENT_TYPE':self.headers['Content-Type']})

        print form
        return
#         # Get relevant POST parameters.
#         # Check parameters.
#         if any(s not in form for s in ('username', 'password', 'nodes')):
#             self.send_error(400, 'Invalid parameters (username, password, nodes).')
#             return
# 
#         try:
#             username = form['username'].value
#             password = form['password'].value
#             num_nodes = int(form['nodes'].value)
#         except:
#             self.send_error(400, 'Invalid parameter values.')
#             return
# 
#         # Check username and password
#         if password != PASSWORD:
#             self.send_error(400, 'Invalid password.%s %s' % (password, PASSWORD))
#             return
# 
#         # Check number of nodes requested
#         if (num_nodes < 1):
#             self.send_error(400, 'Invalid number of requested nodes.')
#             return
# 
#         # Assign IDs to this job.
#         job = 'j' + str(uuid.uuid4())
#         filename = '/home/.maxwell.' + job
#         project = username.replace('@', '-at-')
# 
#         self.log_message("[%s] %s requests %s nodes" % \
#                             (job, username, num_nodes))
# 
#         # Save input data to a file.
#         infile = open(filename + '.i', 'w')
#         shutil.copyfileobj(form['in'].file, infile)
#         infile.close()
# 
#         # Create project.
#         project_file = open(filename + '.' + project, 'w')
#         project_file.write(project_template % (project, 100))
#         project_file.close()
#         try:
#             subprocess.check_output(shlex.split( \
#                         "qconf -Aprj %s" % (filename + '.' + project)), \
# 			stderr=subprocess.STDOUT)
#         except: pass
# 
#         os.remove(filename + '.' + project)
#     
#         # Submit job.
#         qsub_cmd =  "qsub " + \
#                     "-P \"" + project + "\" " \
#                     "-pe orte " + str(num_nodes)  + " " + \
#                     "-N " + "\"" + job + "\"" + " " + \
#                     "-o /home/fdfd_out -e /home/fdfd_err " + \
#                     "-b y mpirun /root/maxwell_env/bin/python " + \
#                         "/home/fdfd/fdfd.py \"" +  filename + "\""
# 
#         # print qsub_cmd 
#         try:
#             response = subprocess.check_output(shlex.split(qsub_cmd))
#             self.log_message("[%s] %s" % (job, response.strip()))
#         except subprocess.CalledProcessError as e:
#             if e.returncode == 1: # Not actually error, just no nodes yet.
#                 pass
#             else:
#                 self.send_error(400, 'Internal server error, job not submitted.')
#                 return
#         except:
#             self.send_error(400, 'Fatal internal error, job not submitted.')
#             return
# 
#         # Header information for user to obtain results.
#         self.send_response(200) # Simple header.
#         self.send_header('Maxwell-Redirect', 'http://maxwell-cloudfront.lightlabs.co') 
#         self.send_header('Maxwell-Name', job)
#         self.end_headers()
# 
#         # Stream residual information until job is done.
#         job_started = lambda: os.path.isfile('/home/' + job + '.start')
#         job_finished = lambda: os.path.isfile('/home/' + job + '.stop')
# 
#         check_state = lambda state: job in \
#             subprocess.check_output(shlex.split("qstat -xml -s " + state))
# 
#         status_file = open(filename + '.s', 'w')
#         status_file = open(filename + '.s', 'r')
#         start_time = time.time() # So we know how long we're waiting.
#         start_logged = False
#         while True:
#             # Put this in front to make sure we always get everything
#             # from the status file.
#             is_done = job_finished() 
# 
#             time.sleep(PAUSE_TIME)
#             if job_started():
#                 if not start_logged:
#                     self.log_message("[%s] started" % (job))
#                     start_logged = True
#                 self.wfile.write("".join(status_file.readlines()))
#             else:
#                 self.wfile.write("WAIT %1.1f seconds\n" % (time.time() - start_time))
# 
#             if is_done:
#                 self.log_message("[%s] finished" % job)
#                 break
# 
#         # Delete input and status files.
#         os.remove(filename + '.i')
#         os.remove(filename + '.s')

    def do_GET(self):
        """ Serves simulation output files.
        
            Basically acts like a file downloader.
            However, to prevent users from downloading whatever they way,
            we check to make sure the file ending corresponds to a valid 
            electromagnetic field.
        """

        filename = "/home" + self.path

        # Check for valid file ending and existence.
        if not((filename[-5:] in ('_real', '_imag')) and \
                (filename[-6] in ('x', 'y', 'z')) and \
                (filename[-7] in ('E', 'H')) and \
                os.path.isfile(filename)):
            self.send_error(404, "file not found")

        else: # File is valid and exists, send it to the user.
            # Rename to prevent re-downloading.
            os.rename(filename, filename + ".send") 

            f = open(filename + ".send", 'r')
            self.send_response(200)
            self.send_header('Content-type', 'text/plain') 
            self.end_headers()
            shutil.copyfileobj(f, self.wfile)
            f.close()

            os.remove(filename + ".send") # Delete file after transmission.


class ForkingHTTPServer(ForkingMixIn, BaseHTTPServer.HTTPServer):
    """ We use a multi-process version of BaseHTTPServer. """


if __name__ == '__main__':
    print "Clearing all jobs out of scheduler..."
    try:
        subprocess.call(shlex.split('qdel -u "*"')) 
    except: pass
    server_address = ("", int(sys.argv[1]))
    print "Serving at", server_address
    print "Password: ", PASSWORD
    httpd = ForkingHTTPServer(server_address, MaxwellHandler)
    httpd.serve_forever()

